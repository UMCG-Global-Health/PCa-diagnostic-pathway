# =============================================================================
# Fast MRI Prostate Cancer Screening Model
# =============================================================================
# Description: Discrete-event simulation (DES) of a prostate cancer (PCa)
#              diagnostic pathway in the Dutch healthcare system, evaluating
#              the health economic impact of AI-assisted MRI reading.
#
# Model structure:
#   1. GP visit + PSA/DRE blood draw
#   2. Urologist intake
#   3. MRI scan + radiologist assessment
#   4. Biopsy + pathologist assessment
#   5. Treatment or active surveillance
#
# Patient risk categories:
#   1 = Clinically significant PCa (Gleason >= 7)
#   2 = Non-significant PCa
#   3 = No prostate cancer
#
# Author:      Pim WM van Dorst
# Affiliation: University Medical Center Groningen
# Date:        20-05-2026
# License:     MIT
# =============================================================================


# 0. SETUP --------------------------------------------------------------------

rm(list = ls())

required_packages <- c("here", "tidyverse", "simmer", "simmer.plot", "ggpubr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

library(here)
library(tidyverse)
library(simmer)
library(simmer.plot)
library(ggpubr)

options(scipen = 10)  # Suppress scientific notation in output


# 1. PARAMETERS ---------------------------------------------------------------

## 1.1 Discount rates ---------------------------------------------------------

discount_rate_qaly <- 1.5  # % per year
discount_rate_cost <- 3.0  # % per year
days_in_year       <- 365.25

# Derived daily rates
daily_rate_qaly <- discount_rate_qaly / days_in_year
daily_rate_cost <- discount_rate_cost / days_in_year


## 1.2 Patient population -----------------------------------------------------

# Probability of each risk group (used for sampling):
#   Prevalence PCa in men: 11.1%
#   Proportion clinically significant (Gleason >= 7): 42%
#   Proportion non-significant: 58%
prob_risk <- c(
  sig    = 0.111 * 0.42,   # Risk group 1: significant PCa
  no_sig = 0.111 * 0.58,   # Risk group 2: non-significant PCa
  no_dis = 1 - 0.111        # Risk group 3: no disease
)

# Pre-sample risk for all patients (vectorised for performance)
n_patients_total <- 277042
set.seed(123)
assign_risk_v <- sample(
  x       = c(1, 2, 3),
  size    = n_patients_total,
  replace = TRUE,
  prob    = prob_risk
)

# Fraction of patients entering the simulation (based on diagnostic pathway)
perc_in <- 0.99 * 0.9 * 0.97 * 0.9 * 0.111

## 1.3 Test characteristics ---------------------------------------------------

# Each test is defined by sensitivity and specificity per risk group.
# Groups: sig (1), no_sig (2), no_dis (3).
# For no_dis, sensitivity is combined across sig and no_sig using the
# addition rule for non-mutually-exclusive events.

### PSA + DRE
psa_sens_cs  <- 0.97
dre_sens_cs  <- 0.286
psa_spec_cs  <- 0.089
dre_spec_cs  <- 0.94
psa_spec_all <- 0.17
dre_spec_all <- 0.94

psa_dre_sens_cs  <- psa_sens_cs + dre_sens_cs - (psa_sens_cs * dre_sens_cs)
psa_dre_spec_cs  <- psa_spec_cs * dre_spec_cs
psa_dre_spec_all <- psa_spec_all * dre_spec_all
psa_dre_sens_any <- 0.46

psa_test_characteristics <- list(
  sig    = list(sensitivity = psa_dre_sens_cs,  specificity = psa_dre_spec_cs),
  no_sig = list(sensitivity = psa_dre_sens_any, specificity = 0.83),
  no_dis = list(
    sensitivity = psa_dre_sens_cs + psa_dre_sens_any - (psa_dre_sens_cs * psa_dre_sens_any),
    specificity = psa_dre_spec_all
  )
)

### Urologist assessment
uro_test_characteristics <- list(
  sig    = list(sensitivity = 0.65, specificity = 0.74),
  no_sig = list(sensitivity = 0.70, specificity = 0.61),
  no_dis = list(
    sensitivity = 0.65 + 0.70 - (0.65 * 0.70),
    specificity = 0.61
  )
)

### MRI (standard)
mri_test_characteristics <- list(
  sig    = list(sensitivity = 0.91, specificity = 0.577),
  no_sig = list(sensitivity = 0.70, specificity = 0.577),
  no_dis = list(
    sensitivity = 0.91 + 0.70 - (0.91 * 0.70),
    specificity = 0.577
  )
)

### MRI (AI-assisted) — applied after year 5 of the simulation
mri_test_characteristics_ai <- list(
  sig    = list(sensitivity = 0.91 * 0.98, specificity = 0.577 * 0.98),
  no_sig = list(sensitivity = 0.70 * 0.98, specificity = 0.577 * 0.98),
  no_dis = list(
    sensitivity = (0.91 * 0.98) + (0.70 * 0.98) - ((0.91 * 0.98) * (0.70 * 0.98)),
    specificity = 0.577 * 0.98
  )
)

### Biopsy
biop_test_characteristics <- list(
  sig    = list(sensitivity = 0.80, specificity = 0.94),
  no_sig = list(sensitivity = 0.51, specificity = 1.00),
  no_dis = list(
    sensitivity = 0.80 + 0.51 - (0.80 * 0.51),
    specificity = 1.00
  )
)


## 1.4 Utility weights --------------------------------------------------------

u_diag_phase <- 0.90  # During diagnostic workup
u_diag       <- 0.80  # When diagnosed
u_trt_prostatect <- (2 * 0.67 + 10 * 0.77) / 12  # Prostatectomy
u_trt_radi       <- (2 * 0.73 + 10 * 0.78) / 12  # Radiotherapy
u_as         <- 0.97  # Active surveillance
u_aftertrt   <- 0.95  # One year post-treatment
u_out        <- 1.00  # Undetected (background health)


## 1.5 Unit costs (euros) -----------------------------------------------------

c_gp_first    <- 33.10
c_blooddraw   <- 31.72
c_urologyconsult <- 145.00
c_mri         <- 283.00
c_biopsy      <- 481 * 1.017 * 1.026 * 1.013 * 1.027 * 1.1 * 1.038 * 1.033  # Inflated
c_treatment   <- 7253 * 1.013 * 1.027 * 1.1 * 1.038 * 1.033                  # Inflated
c_diagnostics <- 1415.00  # Total cost of MRI + biopsy (used in sensitivity)


## 1.6 Resource capacities and service times ----------------------------------

# GP
capacity_gp   <- 1    # Average PCa-relevant consults per working day
t_gp_min      <- 15   # Average  minutes per consult

# Blood draw
capacity_blooddraw <- 10  # Average draws per working day
t_blooddraw_min    <- 5   # Average minutes per draw

# Urologist
capacity_urologist   <- 6   # Average PCa appointments per day
t_urologyconsult_min <- 15  # Average minutes per consult

# MRI scanner
capacity_mri <- 4.4   # Average PCa patients per day per scanner
t_mri_min    <- 25    # Average minutes per scan (+ 5 min prep hardcoded)

# Radiologist
cap_radiologist             <- 4.8  # Average MRI reports per day per radiologist
t_mri_radiologist_case_min  <- 35   # Average minutes per report (normtijden)

# Biopsy (urologist performing)
capacity_biopsy <- 2   # Average biopsies per day

# Pathologist
cap_pathologist        <- 3   # Average biopsies analysed per day
t_biopsy_pathologist_min <- 60  # Average minutes per case

# Waiting times
t_days_psa_result <- 1        # Days until PSA result is returned
t_biopsy_analysis <- 1        # Days before pathologist starts analysis
t_checkbiop       <- 2 * days_in_year  # 2-year follow-up interval after negative biopsy

# Capacity correction factor (working week / calendar week)
cap_recalc <- 7 / 5

# AI efficiency reduction factors (applied after year 5; set to 1 = no change)
# - Strategy 1 - AI-MRI strategy : mri_ai_reduction = 3, rad_ai_reduction = 1
# - Strategy 2 - AI-assessment : mri_ai_reduction = 1, rad_ai_reduction = 2
# - Strategy 3 - Combined AI-MRI & AI-assessment : mri_ai_reduction = 3, rad_ai_reduction = 2
mri_ai_reduction <- 1
rad_ai_reduction <- 1


# 2. HELPER FUNCTIONS ---------------------------------------------------------

## 2.1 Risk assignment --------------------------------------------------------

#' Assign pre-sampled risk group to patient
#'
#' @param number Patient index (extracted from simmer patient name)
#' @param risk_vector Pre-sampled risk vector (default: assign_risk_v)
#' @return Integer: 1 (sig PCa), 2 (non-sig PCa), 3 (no disease)
assign_risk <- function(number, risk_vector = assign_risk_v) {
  risk_vector[[number + 1]]
}


## 2.2 Risk progression -------------------------------------------------------

#' Reassign risk after negative biopsy to account for disease progression
#'
#' Progression probability modelled as exponential with 32% at 15 years
#' (Tosoian et al.).
#'
#' @param risk   Current risk group (1, 2, or 3)
#' @param t_start Start time (days)
#' @param t_now   Current simulation time (days)
#' @return Updated risk group (integer)
reassign_risk <- function(risk, t_start, t_now) {
  years      <- (t_now - t_start) / days_in_year
  progression <- 1 - exp(-(-log(1 - 0.32) / 15) * years)

  if (risk == 2) {
    sample(c(1, 2), size = 1, prob = c(progression, 1 - progression))
  } else {
    risk
  }
}


## 2.3 Mortality functions ----------------------------------------------------

#' Sample survival outcome for untreated/background mortality
#'
#' Returns 1 (alive → continue) or 2 (dead → exit).
#' Significant PCa: 60% 15-year mortality (Albertsen et al.)
#' Background: based on Dutch men aged 45+ (~0.134%/year)
#'
#' @param risk    Risk group (1, 2, 3)
#' @param t_start Start time (days)
#' @param t_now   Current simulation time (days)
#' @return 1 (alive) or 2 (dead)
fun_pca_death <- function(risk, t_start, t_now) {
  years <- (t_now - t_start) / days_in_year

  if (risk == 1) {
    prob_death <- 1 - exp(-(-log(1 - 0.60) / 15) * years)
  } else {
    prob_death <- 1 - exp(-(-log(1 - 0.00134)) * years)
  }

  death <- sample(c(TRUE, FALSE), size = 1, prob = c(prob_death, 1 - prob_death))
  if (death) 2 else 1
}


#' Sample survival outcome during active treatment phase
#'
#' Returns 1 (dead → flag PCa death) or 2 (alive → continue to treatment).
#'
#' @param risk    Risk group (1, 2, 3)
#' @param t_start Start time (days)
#' @param t_now   Current simulation time (days)
#' @return 1 (PCa death) or 2 (alive)
fun_pca_death_pca <- function(risk, t_start, t_now) {
  years <- (t_now - t_start) / days_in_year

  if (risk == 1) {
    prob_death <- 1 - exp(-(-log(1 - 0.60) / 15) * years)
  } else {
    prob_death <- 1 - exp(-(-log(1 - 0.00134)) * years)
  }

  death <- sample(c(TRUE, FALSE), size = 1, prob = c(prob_death, 1 - prob_death))
  if (death) 1 else 2
}


## 2.4 Resource service time functions ----------------------------------------

#' Compute waiting time contributed by GP consults
#' @return Days of waiting time per patient
t_gp <- function(capacity = capacity_gp, t_min = t_gp_min) {
  t_standard  <- 15
  cap_corr    <- t_standard / t_min
  (1 / (capacity * cap_corr)) * cap_recalc
}

#' Compute waiting time contributed by blood draw
#' @return Days of waiting time per patient
t_blooddraw <- function(capacity = capacity_blooddraw, t_min = t_blooddraw_min) {
  t_standard <- 5
  cap_corr   <- t_standard / t_min
  (1 / (capacity * cap_corr)) * cap_recalc
}

#' Compute waiting time contributed by urologist consult
#' @return Days of waiting time per patient
t_urologyconsult <- function(capacity = capacity_urologist, t_min = t_urologyconsult_min) {
  t_standard <- 15
  cap_corr   <- t_standard / t_min
  (1 / (capacity * cap_corr)) * cap_recalc
}

#' Compute waiting time contributed by MRI scan
#' @param t_scan_min Actual scan duration in minutes (may differ from baseline)
#' @return Days of waiting time per patient
t_mri_scan <- function(capacity = capacity_mri, t_scan_min = t_mri_min) {
  t_standard <- t_mri_min + 5  # Scan + preparation
  cap_corr   <- t_standard / (t_scan_min + 5)
  (1 / (capacity * cap_corr)) * cap_recalc
}

#' Compute waiting time contributed by radiologist MRI assessment
#' @param t_case_min Actual radiologist time per case in minutes
#' @return Days of waiting time per patient
t_mri_radiologist_analysis <- function(capacity = cap_radiologist,
                                        t_case_min = t_mri_radiologist_case_min) {
  t_standard <- 35
  cap_corr   <- t_standard / t_case_min
  (1 / (capacity * cap_corr)) * cap_recalc
}

#' Compute waiting time contributed by biopsy
#' @return Days of waiting time per patient
t_biopsy_wait <- function(capacity = capacity_biopsy) {
  (1 / capacity) * cap_recalc
}

#' Sample waiting time before pathologist starts analysis
#' @return Days of waiting time
t_biopsyanalysis_timeout <- function(t_days = t_biopsy_analysis) {
  1 / t_days
}

#' Compute waiting time contributed by pathologist biopsy analysis
#' @return Days of waiting time per patient
t_biopsypathologist <- function(capacity = cap_pathologist,
                                 t_case_min = t_biopsy_pathologist_min) {
  t_standard <- 60
  cap_corr   <- t_standard / t_case_min
  (1 / (capacity * cap_corr)) * cap_recalc
}


## 2.5 Health economic outcome functions --------------------------------------

#' Calculate discounted QALYs over a time interval
#'
#' QALYs are accumulated only within the analysis window (years 5–15 of
#' the simulation). Earlier periods are excluded as burn-in.
#'
#' @param t_start  Start of the health state (simulation days)
#' @param t_now    End of the health state (simulation days)
#' @param u_v      Utility value (0–1)
#' @param disc     Daily discount rate (default: daily_rate_qaly)
#' @return Discounted QALYs (numeric)
fun_qaly <- function(t_start, t_now, u_v, disc = daily_rate_qaly) {
  max_time <- 15 * days_in_year
  burn_in  <- 5  * days_in_year

  t_start_adj <- max(min(t_start, max_time) - burn_in, 0)
  t_now_adj   <- max(min(t_now,   max_time) - burn_in, 0)

  if (t_start >= max_time) return(0)

  time_seq     <- seq(from = t_start_adj, to = t_now_adj, by = 1)
  disc_factors <- 1 / ((1 + disc / 100) ^ time_seq)
  sum((u_v / days_in_year) * disc_factors)
}


#' Calculate discounted QALYs for the treatment phase
#'
#' Applies a three-phase utility profile:
#'   - First 2 months: u = 0.70 (acute treatment)
#'   - Months 2–12:    u = 0.78 (recovery)
#'   - Year 1+:        u = 0.95 (post-treatment)
#'
#' @param t_start Start of treatment (simulation days)
#' @param disc    Daily discount rate (default: daily_rate_qaly)
#' @return Discounted QALYs (numeric)
fun_qaly_trt <- function(t_start, disc = daily_rate_qaly) {
  max_time <- 15 * days_in_year
  burn_in  <- 5  * days_in_year

  cur_time  <- max(min(t_start, max_time) - burn_in, 0)
  diff_time <- (max_time - burn_in) - cur_time

  u_f <- min(diff_time / days_in_year, 2 / 12) * 0.70
  u_s <- if (diff_time / days_in_year <= 2 / 12) {
    0
  } else if (diff_time / days_in_year <= 1) {
    (diff_time / days_in_year - 2 / 12) * 0.78
  } else {
    (10 / 12) * 0.78
  }
  u_t <- if (diff_time / days_in_year > 1) (diff_time / days_in_year - 1) * 0.95 else 0

  u_tot <- u_f + u_s + u_t

  max_now  <- max_time - burn_in
  max_start <- min(cur_time, max_now)
  time_seq  <- seq(from = max_start, to = max_now, by = 1)
  disc_factors <- 1 / ((1 + disc / 100) ^ time_seq)
  sum(u_tot / (max_now - max_start + 1) * disc_factors)
}


#' Calculate discounted cost at a given simulation time point
#'
#' Costs are discounted relative to the start of the analysis period (year 5).
#'
#' @param t_now   Current simulation time (days)
#' @param c_case  Cost amount (euros)
#' @param disc    Daily discount rate (default: daily_rate_cost)
#' @return Discounted cost (euros)
cost_disc <- function(t_now, c_case, disc = daily_rate_cost) {
  max_time <- 15 * days_in_year
  burn_in  <- 5  * days_in_year

  cur_time <- max(min(t_now, max_time) - burn_in, 0)
  time_seq  <- seq(from = 0, to = cur_time, by = 1)
  disc_factors <- 1 / ((1 + disc / 100) ^ time_seq)
  sum(c_case / (cur_time + 1) * disc_factors)
}


# 3. ARRIVAL TIMES ------------------------------------------------------------

# Annual patient volumes based on Dutch cancer registry data.
# First 5 years (2020–2024): observed counts.
# Years 6–16: projected (flat at 2024 level for the no-trend scenario).
# Volumes are divided by 10 for computational tractability.

build_arrivals <- function(volume_per_year, perc_in, scale = 10) {
  volumes <- volume_per_year / perc_in / scale
  times   <- vector("list", length(volumes))

  for (i in seq_along(volumes)) {
    start <- if (i == 1) 0 else tail(times[[i - 1]], 1)
    times[[i]] <- seq(start, i * days_in_year, length.out = volumes[[i]])
  }
  unlist(times)
}

# Observed volumes 2020–2024
observed_volumes <- c(12541, 13714, 14698, 14422, 15212) * 0.955

# Projected volumes 2025–2035 (two scenarios)
projected_trend  <- c(14211, 14492, 14773, 15054, 15335, 15617, 15898, 16179, 16460, 16741, 17022)
projected_flat   <- rep(15212 * 0.955, 11)

arrival_times_trend <- build_arrivals(
  volume_per_year = c(observed_volumes, projected_trend),
  perc_in = perc_in
)

arrival_times_flat <- build_arrivals(
  volume_per_year = c(observed_volumes, projected_flat),
  perc_in = perc_in
)


# 4. PATIENT TRAJECTORY -------------------------------------------------------

traj_treat_new <- trajectory(name = "prostate_cancer_pathway") %>%

  # --- Initialise patient attributes ---
  set_attribute("initialize",      1) %>%
  set_attribute("initialize_time", function() now(.env = sim_new)) %>%
  set_attribute("risk", function() {
    assign_risk(as.integer(gsub("\\D", "", get_name(.env = sim_new))))
  }) %>%

  # =========================================================================
  # STEP 1: GP VISIT
  # =========================================================================
  set_attribute("queue_GP", function() get_queue_count(.env = sim_new, "GP")) %>%
  seize("GP", 1) %>%
  set_attribute("Costs_gp", function() cost_disc(now(.env = sim_new), c_gp_first)) %>%
  timeout(function() t_gp()) %>%
  release("GP", 1) %>%

  # =========================================================================
  # STEP 2: BLOOD DRAW + PSA RESULT
  # =========================================================================
  set_attribute("queue_BD", function() get_queue_count(.env = sim_new, "Blooddraw"),
                tag = "blooddraw") %>%
  set_attribute("start_bd_time", function() now(.env = sim_new)) %>%
  seize("Blooddraw", 1) %>%
  set_attribute("Costs_bd", function() cost_disc(now(.env = sim_new), c_blooddraw)) %>%
  set_attribute("Start_blooddraw", function() {
    n_bd  <- get_attribute(.env = sim_new, "Start_blooddraw")
    n_bd  <- ifelse(is.na(n_bd), 0, n_bd)
    n_uro <- get_attribute(.env = sim_new, "Start_urologist")
    n_uro <- ifelse(is.na(n_uro), 0, n_uro)
    ifelse(n_bd <= n_uro, n_uro + 1, n_bd + 1)
  }) %>%
  timeout(function() t_blooddraw()) %>%
  release("Blooddraw") %>%
  timeout(t_days_psa_result) %>%

  # =========================================================================
  # STEP 3: PSA / DRE RESULT — branch to urologist or exit
  # =========================================================================
  branch(
    function() {
      risk      <- get_attribute(.env = sim_new, "risk")
      sens      <- psa_test_characteristics[[risk]]$sensitivity
      spec      <- psa_test_characteristics[[risk]]$specificity
      has_dis   <- risk %in% c(1, 2)
      result    <- sample(c(TRUE, FALSE), 1,
                          prob = if (has_dis) c(sens, 1 - sens) else c(1 - spec, spec))
      if (result) 2 else 1
    },
    continue = c(FALSE, FALSE),

    # --- Branch A: Negative PSA — monitor for up to 10 years ---
    trajectory("no_treatment_after_psa") %>%
      set_attribute("no_trt_after_PSA_done", 1, mod = "+") %>%
      set_attribute("Stop", 1, mod = "+") %>%
      set_attribute("QALY", function() {
        fun_qaly(get_attribute(.env = sim_new, "start_bd_time"), now(.env = sim_new), u_diag_phase)
      }) %>%
      set_attribute("no_trt_loop_psa", 1, mod = "+", tag = "no_trt_loop_psa") %>%
      set_attribute("no_trt_start", function() {
        n <- get_attribute(.env = sim_new, "no_trt_loop_psa")
        n <- ifelse(is.na(n), 0, n)
        if (n <= 1) get_attribute(.env = sim_new, "start_bd_time") else now(.env = sim_new)
      }) %>%
      timeout(1 * days_in_year) %>%
      branch(
        option = function() fun_pca_death(
          get_attribute(.env = sim_new, "risk"),
          get_attribute(.env = sim_new, "no_trt_start"),
          now(.env = sim_new)
        ),
        continue = c(FALSE, FALSE),
        trajectory("alive_no_psa") %>%
          set_attribute("risk", function() reassign_risk(
            get_attribute(.env = sim_new, "risk"),
            get_attribute(.env = sim_new, "no_trt_start"),
            now(.env = sim_new)
          )) %>%
          set_attribute("QALY", function() {
            fun_qaly(get_attribute(.env = sim_new, "no_trt_start"), now(.env = sim_new), u_out)
          }) %>%
          rollback("no_trt_loop_psa", times = 9),
        trajectory("pca_death_no_psa") %>%
          set_attribute("QALY", function() {
            fun_qaly(get_attribute(.env = sim_new, "no_trt_start"),
                     now(.env = sim_new) - 0.5 * days_in_year, u_out)
          }) %>%
          set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_blooddraw"))
      ),

    # --- Branch B: Positive PSA — refer to urologist ---
    trajectory("to_urologist") %>%

      # =====================================================================
      # STEP 4: UROLOGIST INTAKE
      # =====================================================================
      set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist"),
                    tag = "urologist_intake") %>%
      set_attribute("start_from_urol", function() {
        n_start <- get_attribute(.env = sim_new, "Start_blooddraw")
        n_stop  <- get_attribute(.env = sim_new, "Stop")
        n_stop  <- ifelse(is.na(n_stop), 0, n_stop)
        if (n_stop >= n_start) now(.env = sim_new) else get_attribute(.env = sim_new, "start_bd_time")
      }) %>%
      set_attribute("Start_urologist", function() {
        n <- get_attribute(.env = sim_new, "Start_urologist")
        n <- ifelse(is.na(n), 0, n)
        n + 1
      }) %>%
      seize("Urologist", 1) %>%
      set_attribute("Cost_urologist", function() cost_disc(now(.env = sim_new), c_urologyconsult)) %>%
      timeout(function() t_urologyconsult()) %>%
      release("Urologist", 1) %>%

      # =====================================================================
      # STEP 5: UROLOGIST ASSESSMENT — branch to MRI or exit
      # =====================================================================
      branch(
        function() {
          risk    <- get_attribute(.env = sim_new, "risk")
          sens    <- uro_test_characteristics[[risk]]$sensitivity
          spec    <- uro_test_characteristics[[risk]]$specificity
          has_dis <- risk %in% c(1, 2)
          result  <- sample(c(TRUE, FALSE), 1,
                            prob = if (has_dis) c(sens, 1 - sens) else c(1 - spec, spec))
          if (result) 2 else 1
        },
        continue = c(FALSE, FALSE),

        # --- Branch A: Negative urologist — monitor ---
        trajectory("no_treatment_after_urologist") %>%
          set_attribute("no_trt_after_urologist_intake_done", 1, mod = "+") %>%
          set_attribute("Stop", 1, mod = "+") %>%
          set_attribute("QALY", function() {
            fun_qaly(get_attribute(.env = sim_new, "start_from_urol"), now(.env = sim_new), u_diag_phase)
          }) %>%
          branch(
            function() {
              n <- get_attribute(.env = sim_new, "no_trt_after_urologist_intake_done")
              if (n == 1) 1 else 2
            },
            continue = c(FALSE, FALSE),
            # First exit: 2-year recheck then possible rollback to blood draw
            trajectory("reassessment_uro") %>%
              set_attribute("start_wait", function() now(.env = sim_new)) %>%
              timeout(2 * days_in_year) %>%
              branch(
                option = function() fun_pca_death(
                  get_attribute(.env = sim_new, "risk"),
                  get_attribute(.env = sim_new, "start_from_urol"),
                  now(.env = sim_new)
                ),
                continue = c(FALSE, FALSE),
                trajectory("re_evaluate_uro") %>%
                  set_attribute("risk", function() reassign_risk(
                    get_attribute(.env = sim_new, "risk"),
                    get_attribute(.env = sim_new, "start_from_urol"),
                    now(.env = sim_new)
                  )) %>%
                  set_attribute("QALY", function() {
                    fun_qaly(get_attribute(.env = sim_new, "start_wait"), now(.env = sim_new), u_out)
                  }) %>%
                  rollback("blooddraw", times = 1),
                trajectory("pca_death_uro") %>%
                  set_attribute("QALY", function() {
                    fun_qaly(get_attribute(.env = sim_new, "start_wait"),
                             now(.env = sim_new) - days_in_year, u_out)
                  }) %>%
                  set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
              ),
            # Subsequent exits: annual monitoring loop (up to 9 iterations)
            trajectory("monitoring_loop_uro") %>%
              set_attribute("no_trt_loop_uro", 1, mod = "+", tag = "no_trt_loop_uro") %>%
              set_attribute("no_trt_start", function() now(.env = sim_new)) %>%
              set_attribute("add_time_start", function() {
                n <- get_attribute(.env = sim_new, "no_trt_loop_uro")
                if (n == 1) get_attribute(.env = sim_new, "start_from_urol") else now(.env = sim_new)
              }) %>%
              timeout(1 * days_in_year) %>%
              branch(
                option = function() fun_pca_death(
                  get_attribute(.env = sim_new, "risk"),
                  get_attribute(.env = sim_new, "add_time_start"),
                  now(.env = sim_new)
                ),
                continue = c(FALSE, FALSE),
                trajectory("alive_monitoring_uro") %>%
                  set_attribute("risk", function() reassign_risk(
                    get_attribute(.env = sim_new, "risk"),
                    get_attribute(.env = sim_new, "add_time_start"),
                    now(.env = sim_new)
                  )) %>%
                  set_attribute("QALY", function() {
                    fun_qaly(get_attribute(.env = sim_new, "no_trt_start"), now(.env = sim_new), u_out)
                  }) %>%
                  rollback("no_trt_loop_uro", times = 9),
                trajectory("pca_death_monitoring_uro") %>%
                  set_attribute("QALY", function() {
                    fun_qaly(get_attribute(.env = sim_new, "no_trt_start"),
                             now(.env = sim_new) - 0.5 * days_in_year, u_out)
                  }) %>%
                  set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
              )
          ),

        # --- Branch B: Positive urologist — proceed to MRI ---
        trajectory("to_mri") %>%

          # ==================================================================
          # STEP 6: MRI SCAN + RADIOLOGIST ASSESSMENT
          # ==================================================================
          set_attribute("queue_MRI", function() get_queue_count(.env = sim_new, "MRI")) %>%
          seize("MRI", 1) %>%
          set_attribute("Start_mri", function() get_attribute(.env = sim_new, "Start_urologist")) %>%
          set_attribute("Cost_mri", function() cost_disc(now(.env = sim_new), c_mri)) %>%
          timeout(function() {
            t_scan <- if (now(.env = sim_new) <= 5 * days_in_year) t_mri_min else t_mri_min / mri_ai_reduction
            t_mri_scan(t_scan_min = t_scan)
          }) %>%
          release("MRI", 1) %>%
          set_attribute("queue_Radiologist", function() get_queue_count(.env = sim_new, "Radiologist")) %>%
          seize("Radiologist", 1) %>%
          timeout(function() {
            t_rad <- if (now(.env = sim_new) <= 5 * days_in_year) {
              t_mri_radiologist_case_min
            } else {
              t_mri_radiologist_case_min / rad_ai_reduction
            }
            t_mri_radiologist_analysis(t_case_min = t_rad)
          }) %>%
          release("Radiologist", 1) %>%

          # Urologist consult to discuss MRI results
          set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
          seize("Urologist", 1) %>%
          set_attribute("Cost_urologist", function() cost_disc(now(.env = sim_new), c_urologyconsult)) %>%
          timeout(function() t_urologyconsult()) %>%
          release("Urologist") %>%

          # ==================================================================
          # STEP 7: MRI RESULT — branch to biopsy or exit
          # ==================================================================
          branch(
            function() {
              risk    <- get_attribute(.env = sim_new, "risk")
              chars   <- if (now(.env = sim_new) <= 5 * days_in_year) {
                mri_test_characteristics[[risk]]
              } else {
                mri_test_characteristics_ai[[risk]]
              }
              has_dis <- risk %in% c(1, 2)
              result  <- sample(c(TRUE, FALSE), 1,
                                prob = if (has_dis) {
                                  c(chars$sensitivity, 1 - chars$sensitivity)
                                } else {
                                  c(1 - chars$specificity, chars$specificity)
                                })
              if (result) 2 else 1
            },
            continue = c(FALSE, FALSE),

            # --- Branch A: Negative MRI (PI-RADS 1–2) — monitor ---
            trajectory("no_treatment_after_mri") %>%
              set_attribute("no_trt_after_mri_done", 1, mod = "+") %>%
              set_attribute("Stop", 1, mod = "+") %>%
              set_attribute("QALY", function() {
                fun_qaly(get_attribute(.env = sim_new, "start_from_urol"), now(.env = sim_new), u_diag_phase)
              }) %>%
              branch(
                function() {
                  n <- get_attribute(.env = sim_new, "no_trt_after_mri_done")
                  if (n == 1) 1 else 2
                },
                continue = c(FALSE, FALSE),
                # First negative MRI: 2-year recheck
                trajectory("reassessment_mri") %>%
                  set_attribute("start_wait", function() now(.env = sim_new)) %>%
                  timeout(2 * days_in_year) %>%
                  branch(
                    option = function() fun_pca_death(
                      get_attribute(.env = sim_new, "risk"),
                      get_attribute(.env = sim_new, "start_from_urol"),
                      now(.env = sim_new)
                    ),
                    continue = c(FALSE, FALSE),
                    trajectory("re_evaluate_mri") %>%
                      set_attribute("risk", function() reassign_risk(
                        get_attribute(.env = sim_new, "risk"),
                        get_attribute(.env = sim_new, "start_from_urol"),
                        now(.env = sim_new)
                      )) %>%
                      set_attribute("QALY", function() {
                        fun_qaly(get_attribute(.env = sim_new, "start_wait"), now(.env = sim_new), u_out)
                      }) %>%
                      rollback("blooddraw", times = 1),
                    trajectory("pca_death_mri") %>%
                      set_attribute("QALY", function() {
                        fun_qaly(get_attribute(.env = sim_new, "start_wait"),
                                 now(.env = sim_new) - days_in_year, u_out)
                      }) %>%
                      set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
                  ),
                # Subsequent: annual monitoring loop (up to 9 iterations)
                trajectory("monitoring_loop_mri") %>%
                  set_attribute("no_trt_loopmri", 1, mod = "+", tag = "no_trt_loopmri") %>%
                  set_attribute("no_trt_start", function() now(.env = sim_new)) %>%
                  set_attribute("add_time_start", function() {
                    n <- get_attribute(.env = sim_new, "no_trt_loopmri")
                    if (n == 1) get_attribute(.env = sim_new, "start_from_urol") else now(.env = sim_new)
                  }) %>%
                  timeout(1 * days_in_year) %>%
                  branch(
                    option = function() fun_pca_death(
                      get_attribute(.env = sim_new, "risk"),
                      get_attribute(.env = sim_new, "add_time_start"),
                      now(.env = sim_new)
                    ),
                    continue = c(FALSE, FALSE),
                    trajectory("alive_monitoring_mri") %>%
                      set_attribute("risk", function() reassign_risk(
                        get_attribute(.env = sim_new, "risk"),
                        get_attribute(.env = sim_new, "add_time_start"),
                        now(.env = sim_new)
                      )) %>%
                      set_attribute("QALY", function() {
                        fun_qaly(get_attribute(.env = sim_new, "no_trt_start"), now(.env = sim_new), u_out)
                      }) %>%
                      rollback("no_trt_loopmri", times = 9),
                    trajectory("pca_death_monitoring_mri") %>%
                      set_attribute("QALY", function() {
                        fun_qaly(get_attribute(.env = sim_new, "no_trt_start"),
                                 now(.env = sim_new) - 0.5 * days_in_year, u_out)
                      }) %>%
                      set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
                  )
              ),

            # --- Branch B: Positive MRI (PI-RADS >= 3) — proceed to biopsy ---
            trajectory("to_biopsy") %>%
              set_attribute("go_to_biopsy", 1, mod = "+") %>%

              # ================================================================
              # STEP 8: BIOPSY
              # ================================================================
              set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
              seize("Urologist", 1) %>%
              set_attribute("Start_biopsy", function() get_attribute(.env = sim_new, "Start_urologist")) %>%
              set_attribute("Cost_biopsy", function() cost_disc(now(.env = sim_new), c_biopsy)) %>%
              timeout(function() t_urologyconsult(t_min = 2 * t_urologyconsult_min)) %>%
              release("Urologist", 1) %>%
              timeout(function() t_biopsyanalysis_timeout()) %>%
              set_attribute("queue_Pathologist", function() get_queue_count(.env = sim_new, "Pathologist")) %>%
              seize("Pathologist", 1) %>%
              timeout(function() t_biopsypathologist()) %>%
              release("Pathologist", 1) %>%

              # Urologist consult to discuss biopsy result
              set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
              seize("Urologist", 1) %>%
              set_attribute("Cost_urologist", function() cost_disc(now(.env = sim_new), c_urologyconsult)) %>%
              timeout(function() t_urologyconsult()) %>%
              release("Urologist", 1) %>%

              # ==============================================================
              # STEP 9: BIOPSY RESULT — branch to treatment or exit
              # ==============================================================
              branch(
                function() {
                  risk    <- get_attribute(.env = sim_new, "risk")
                  sens    <- biop_test_characteristics[[risk]]$sensitivity
                  spec    <- biop_test_characteristics[[risk]]$specificity
                  has_dis <- risk %in% c(1, 2)
                  result  <- sample(c(TRUE, FALSE), 1,
                                    prob = if (has_dis) c(sens, 1 - sens) else c(1 - spec, spec))
                  if (result) 2 else 1
                },
                continue = c(FALSE, FALSE),

                # --- Branch A: Negative biopsy — 2-year recheck then monitoring ---
                trajectory("no_treatment_after_biopsy") %>%
                  set_attribute("no_trt_after_one_biopsy_done", 1, mod = "+") %>%
                  set_attribute("Stop", 1, mod = "+") %>%
                  set_attribute("QALY", function() {
                    fun_qaly(get_attribute(.env = sim_new, "start_from_urol"), now(.env = sim_new), u_diag_phase)
                  }) %>%
                  branch(
                    function() {
                      n <- get_attribute(.env = sim_new, "no_trt_after_one_biopsy_done")
                      if (n == 1) 1 else 2
                    },
                    continue = c(FALSE, FALSE),
                    # First negative biopsy: 2-year recheck, possible rollback to urologist
                    trajectory("reassessment_biop") %>%
                      set_attribute("no_trt_after_biopt_time", function() now(.env = sim_new)) %>%
                      timeout(t_checkbiop) %>%
                      branch(
                        option = function() fun_pca_death(
                          get_attribute(.env = sim_new, "risk"),
                          get_attribute(.env = sim_new, "start_from_urol"),
                          now(.env = sim_new)
                        ),
                        continue = c(FALSE, FALSE),
                        trajectory("re_evaluate_biop") %>%
                          set_attribute("risk", function() reassign_risk(
                            get_attribute(.env = sim_new, "risk"),
                            get_attribute(.env = sim_new, "start_from_urol"),
                            now(.env = sim_new)
                          )) %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "no_trt_after_biopt_time"),
                                     now(.env = sim_new), u_out)
                          }) %>%
                          rollback("urologist_intake", times = 1),
                        trajectory("pca_death_biop") %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "no_trt_after_biopt_time"),
                                     now(.env = sim_new) - days_in_year, u_out)
                          }) %>%
                          set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
                      ),
                    # Subsequent: annual monitoring loop (up to 8 iterations)
                    trajectory("monitoring_loop_biop") %>%
                      set_attribute("no_trt_loop_biop", 1, mod = "+", tag = "no_trt_loop_biop") %>%
                      set_attribute("no_trt_start", function() now(.env = sim_new)) %>%
                      set_attribute("add_time_start", function() {
                        n <- get_attribute(.env = sim_new, "no_trt_loop_biop")
                        if (n == 1) get_attribute(.env = sim_new, "start_from_urol") else now(.env = sim_new)
                      }) %>%
                      timeout(1 * days_in_year) %>%
                      branch(
                        option = function() fun_pca_death(
                          get_attribute(.env = sim_new, "risk"),
                          get_attribute(.env = sim_new, "add_time_start"),
                          now(.env = sim_new)
                        ),
                        continue = c(FALSE, FALSE),
                        trajectory("alive_monitoring_biop") %>%
                          set_attribute("risk", function() reassign_risk(
                            get_attribute(.env = sim_new, "risk"),
                            get_attribute(.env = sim_new, "add_time_start"),
                            now(.env = sim_new)
                          )) %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "no_trt_start"), now(.env = sim_new), u_out)
                          }) %>%
                          rollback("no_trt_loop_biop", times = 8),
                        trajectory("pca_death_monitoring_biop") %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "no_trt_start"),
                                     now(.env = sim_new) - 0.5 * days_in_year, u_out)
                          }) %>%
                          set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist"))
                      )
                  ),

                # --- Branch B: Positive biopsy — treatment or active surveillance ---
                trajectory("eligible_for_treatment") %>%
                  set_attribute("Treatment_after_MRI_biopsy", 1, mod = "+") %>%

                  # Check for PCa death before treatment starts
                  branch(
                    option = function() fun_pca_death_pca(
                      get_attribute(.env = sim_new, "risk"),
                      get_attribute(.env = sim_new, "start_from_urol"),
                      now(.env = sim_new)
                    ),
                    continue = c(FALSE, FALSE),

                    trajectory("pca_death_before_treatment") %>%
                      set_attribute("QALY", function() {
                        fun_qaly(get_attribute(.env = sim_new, "start_from_urol"),
                                 now(.env = sim_new), u_diag_phase / 2)
                      }) %>%
                      set_attribute("PCa_death", function() get_attribute(.env = sim_new, "Start_urologist")),

                    trajectory("proceed_to_treatment") %>%
                      branch(
                        function() {
                          risk <- get_attribute(.env = sim_new, "risk")
                          if (risk == 1) 1 else 2  # Sig PCa → treatment; non-sig → active surveillance
                        },
                        continue = c(FALSE, FALSE),

                        # --- Curative treatment (significant PCa) ---
                        trajectory("curative_treatment") %>%
                          set_attribute("Start_treatment", function() get_attribute(.env = sim_new, "Start_urologist")) %>%
                          set_attribute("Costs_treatment", function() cost_disc(now(.env = sim_new), c_treatment)) %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "start_from_urol"), now(.env = sim_new), u_diag_phase)
                          }) %>%
                          set_attribute("QALY", function() fun_qaly_trt(t_start = now(.env = sim_new))),

                        # --- Active surveillance (non-significant PCa) ---
                        trajectory("active_surveillance_entry") %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "start_from_urol"), now(.env = sim_new), u_diag_phase)
                          }) %>%
                          set_attribute("Start_as", function() get_attribute(.env = sim_new, "Start_urologist")) %>%
                          set_attribute("starting_as_time", function() now(.env = sim_new)) %>%
                          timeout(0.5 * days_in_year) %>%
                          set_attribute("QALY", function() {
                            fun_qaly(get_attribute(.env = sim_new, "starting_as_time"), now(.env = sim_new), u_as)
                          }) %>%

                          # Active surveillance monitoring loop (up to 8 cycles)
                          set_attribute("active_surveillance", 1, mod = "+", tag = "active_surveillance") %>%
                          set_attribute("back_to_as", 1, mod = "+") %>%
                          set_attribute("start_as_diag_time", function() now(.env = sim_new)) %>%

                          # Blood draw during AS
                          set_attribute("queue_BD", function() get_queue_count(.env = sim_new, "Blooddraw")) %>%
                          seize("Blooddraw", 1) %>%
                          timeout(function() t_blooddraw()) %>%
                          release("Blooddraw") %>%
                          timeout(t_days_psa_result) %>%

                          # Urologist consult with PSA result
                          set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
                          seize("Urologist", 1) %>%
                          timeout(function() t_urologyconsult()) %>%
                          release("Urologist", 1) %>%

                          # PSA/DRE result during AS
                          branch(
                            function() {
                              risk    <- get_attribute(.env = sim_new, "risk")
                              sens    <- psa_test_characteristics[[risk]]$sensitivity
                              spec    <- psa_test_characteristics[[risk]]$specificity
                              has_dis <- risk %in% c(1, 2)
                              result  <- sample(c(TRUE, FALSE), 1,
                                                prob = if (has_dis) c(sens, 1 - sens) else c(1 - spec, spec))
                              # Force periodic MRI at cycles 2, 5, and 7
                              n_as <- get_attribute(.env = sim_new, "back_to_as")
                              result2 <- if (n_as %in% c(2, 5, 7)) TRUE else result
                              if (result2) 2 else 1
                            },
                            continue = c(FALSE, FALSE),

                            # Continue AS without further imaging
                            trajectory("as_continue") %>%
                              set_attribute("Costs_ur1", function() cost_disc(now(.env = sim_new), c_urologyconsult)) %>%
                              set_attribute("QALY", function() {
                                fun_qaly(get_attribute(.env = sim_new, "start_as_diag_time"), now(.env = sim_new), u_diag_phase)
                              }) %>%
                              set_attribute("start_as_wait_time", function() now(.env = sim_new)) %>%
                              set_attribute("time_out_as", function() {
                                n <- get_attribute(.env = sim_new, "back_to_as")
                                ifelse(n < 5, 0.5, 1) * days_in_year
                              }) %>%
                              timeout(function() get_attribute(.env = sim_new, "time_out_as")) %>%
                              set_attribute("QALY", function() {
                                fun_qaly(get_attribute(.env = sim_new, "start_as_wait_time"), now(.env = sim_new), u_as)
                              }) %>%
                              set_attribute("risk", function() reassign_risk(
                                get_attribute(.env = sim_new, "risk"),
                                get_attribute(.env = sim_new, "start_as_diag_time"),
                                now(.env = sim_new)
                              )) %>%
                              rollback("active_surveillance", times = 8),

                            # MRI + biopsy during AS
                            trajectory("as_mri_biopsy") %>%
                              set_attribute("queue_MRI", function() get_queue_count(.env = sim_new, "MRI")) %>%
                              seize("MRI", 1) %>%
                              set_attribute("Costs_diagnostics_as", function() {
                                cost_disc(now(.env = sim_new), c_urologyconsult + c_mri)
                              }) %>%
                              timeout(function() {
                                t_scan <- if (now(.env = sim_new) <= 5 * days_in_year) t_mri_min else t_mri_min / mri_ai_reduction
                                t_mri_scan(t_scan_min = t_scan)
                              }) %>%
                              release("MRI", 1) %>%
                              set_attribute("queue_Radiologist", function() get_queue_count(.env = sim_new, "Radiologist")) %>%
                              seize("Radiologist", 1) %>%
                              timeout(function() {
                                t_rad <- if (now(.env = sim_new) <= 5 * days_in_year) {
                                  t_mri_radiologist_case_min
                                } else {
                                  t_mri_radiologist_case_min / rad_ai_reduction
                                }
                                t_mri_radiologist_analysis(t_case_min = t_rad)
                              }) %>%
                              release("Radiologist", 1) %>%

                              # MRI result during AS
                              branch(
                                function() {
                                  risk  <- get_attribute(.env = sim_new, "risk")
                                  chars <- if (now(.env = sim_new) <= 5 * days_in_year) {
                                    mri_test_characteristics[[risk]]
                                  } else {
                                    mri_test_characteristics_ai[[risk]]
                                  }
                                  has_dis <- risk %in% c(1, 2)
                                  result  <- sample(c(TRUE, FALSE), 1,
                                                    prob = if (has_dis) {
                                                      c(chars$sensitivity, 1 - chars$sensitivity)
                                                    } else {
                                                      c(1 - chars$specificity, chars$specificity)
                                                    })
                                  if (result) 2 else 1
                                },
                                continue = c(FALSE, FALSE),

                                # Negative MRI during AS — return to surveillance
                                trajectory("as_negative_mri") %>%
                                  set_attribute("no_trt_after_mri_done", 1, mod = "+") %>%
                                  set_attribute("Stop_as", 1, mod = "+") %>%
                                  set_attribute("QALY", function() {
                                    fun_qaly(get_attribute(.env = sim_new, "start_as_diag_time"), now(.env = sim_new), u_diag_phase)
                                  }) %>%
                                  set_attribute("start_as_wait_time", function() now(.env = sim_new)) %>%
                                  set_attribute("time_out_as", function() {
                                    n <- get_attribute(.env = sim_new, "back_to_as")
                                    ifelse(n < 5, 0.5, 1) * days_in_year
                                  }) %>%
                                  timeout(function() get_attribute(.env = sim_new, "time_out_as")) %>%
                                  set_attribute("QALY", function() {
                                    fun_qaly(get_attribute(.env = sim_new, "start_as_wait_time"), now(.env = sim_new), u_as)
                                  }) %>%
                                  set_attribute("risk", function() reassign_risk(
                                    get_attribute(.env = sim_new, "risk"),
                                    get_attribute(.env = sim_new, "start_as_diag_time"),
                                    now(.env = sim_new)
                                  )) %>%
                                  rollback("active_surveillance", times = 10),

                                # Positive MRI during AS — biopsy
                                trajectory("as_biopsy") %>%
                                  set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
                                  seize("Urologist", 1) %>%
                                  timeout(function() t_urologyconsult(t_min = 2 * t_urologyconsult_min)) %>%
                                  release("Urologist", 1) %>%
                                  timeout(function() t_biopsyanalysis_timeout()) %>%
                                  set_attribute("queue_Pathologist", function() get_queue_count(.env = sim_new, "Pathologist")) %>%
                                  seize("Pathologist", 1) %>%
                                  timeout(function() t_biopsypathologist()) %>%
                                  release("Pathologist", 1) %>%
                                  set_attribute("queue_Urologist", function() get_queue_count(.env = sim_new, "Urologist")) %>%
                                  seize("Urologist", 1) %>%
                                  timeout(function() t_urologyconsult()) %>%
                                  release("Urologist", 1) %>%
                                  set_attribute("Costs_diagnostics_as", function() cost_disc(now(.env = sim_new), c_biopsy)) %>%

                                  # Biopsy result during AS
                                  branch(
                                    function() {
                                      risk    <- get_attribute(.env = sim_new, "risk")
                                      sens    <- biop_test_characteristics[[risk]]$sensitivity
                                      spec    <- biop_test_characteristics[[risk]]$specificity
                                      has_dis <- risk %in% c(1, 2)
                                      result  <- sample(c(TRUE, FALSE), 1,
                                                        prob = if (has_dis) c(sens, 1 - sens) else c(1 - spec, spec))
                                      if (result) 2 else 1
                                    },
                                    continue = c(FALSE, FALSE),

                                    # Back to AS after negative biopsy (up to 6 cycles)
                                    trajectory("as_back_to_as") %>%
                                      set_attribute("QALY", function() {
                                        fun_qaly(get_attribute(.env = sim_new, "start_as_diag_time"), now(.env = sim_new), u_diag_phase)
                                      }) %>%
                                      set_attribute("time_out_as", function() {
                                        n <- get_attribute(.env = sim_new, "back_to_as")
                                        ifelse(n < 5, 0.5, 1) * days_in_year
                                      }) %>%
                                      set_attribute("start_as_wait_time", function() now(.env = sim_new)) %>%
                                      timeout(function() get_attribute(.env = sim_new, "time_out_as")) %>%
                                      set_attribute("QALY", function() {
                                        fun_qaly(get_attribute(.env = sim_new, "start_as_wait_time"), now(.env = sim_new), u_as)
                                      }) %>%
                                      set_attribute("risk", function() reassign_risk(
                                        get_attribute(.env = sim_new, "risk"),
                                        get_attribute(.env = sim_new, "start_as_diag_time"),
                                        now(.env = sim_new)
                                      )) %>%
                                      rollback("active_surveillance", times = 6),

                                    # Positive biopsy during AS — escalate to treatment
                                    trajectory("as_escalate_to_treatment") %>%
                                      set_attribute("Start_treatment_as", 1, mod = "+") %>%
                                      set_attribute("QALY", function() {
                                        fun_qaly(get_attribute(.env = sim_new, "start_as_diag_time"), now(.env = sim_new), u_diag_phase)
                                      }) %>%
                                      set_attribute("QALY", function() fun_qaly_trt(t_start = now(.env = sim_new))) %>%
                                      set_attribute("Costs_treatment", function() cost_disc(now(.env = sim_new), c_treatment))
                                  )
                              )
                          )
                      )
                  )
              )
          )
      )
  )


# 5. RUN SIMULATION -----------------------------------------------------------

set.seed(123)

sim_new <- simmer() %>%
  add_resource("GP",          capacity = 80)  %>%
  add_resource("Blooddraw",   capacity = 15)  %>%
  add_resource("Urologist",   capacity = 44)  %>%
  add_resource("Radiologist", capacity = 11)  %>%
  add_resource("Pathologist", capacity = 9)   %>%
  add_resource("MRI",         capacity = 12)  %>%
  add_generator(
    name_prefix  = "patient",
    trajectory   = traj_treat_new,
    distribution = at(arrival_times_trend),
    mon          = 2
  )

sim_new %>% reset() %>% run()


# 6. POST-PROCESSING ----------------------------------------------------------

df_attributes <- get_mon_attributes(sim_new)


## 6.1 Queue sizes over time --------------------------------------------------

df_queue <- df_attributes %>%
  filter(grepl("queue", key)) %>%
  dplyr::select(-name, -replication) %>%
  filter(time <= days_in_year * 15) %>%
  group_by(time, key) %>%
  arrange(time) %>%
  slice(n()) %>%
  mutate(key = recode(key,
    "queue_BD"          = "Queue for Blood Draw",
    "queue_GP"          = "Queue for GP",
    "queue_MRI"         = "Queue for MRI",
    "queue_Pathologist" = "Queue for Pathologist",
    "queue_Radiologist" = "Queue for Radiologist",
    "queue_Urologist"   = "Queue for Urologist"
  ))


## 6.2 Diagnostic pathway duration --------------------------------------------

# Exclude burn-in patients (first 5 simulation years)
df_attributes_filtered <- df_attributes %>%
  group_by(name) %>%
  filter(!any(time <= 5 * days_in_year)) %>%
  ungroup()

df_time_diag <- df_attributes_filtered %>%
  filter(time > 5 * days_in_year & time <= days_in_year * 16) %>%
  filter(grepl("Start_treatment|Start_blooddraw|Start_urologist|Stop|PCa_death|Start_as", key)) %>%
  mutate(time = ceiling(time)) %>%
  group_by(value) %>%
  pivot_wider(names_from = "key", values_from = "time") %>%
  mutate(
    start_diag = coalesce(Start_blooddraw, Start_urologist),
    start_trt  = coalesce(Start_as, Start_treatment),
    Stop       = ifelse(is.na(Stop) & !is.na(PCa_death), PCa_death, Stop),
    end_diag   = case_when(
      is.na(Stop) & is.na(start_trt) ~ 15 * days_in_year,
      is.na(Stop)                     ~ start_trt,
      TRUE                            ~ Stop
    ),
    dur_diag = end_diag - start_diag
  )

# Year when diagnostic waiting time first exceeds 21 days
patients_exceeding_21d <- df_time_diag %>% filter(dur_diag > 21)
year_threshold_crossed  <- (min(patients_exceeding_21d$start_diag) - 5 * days_in_year) / days_in_year
cat("Year (relative to analysis start) when 21-day threshold is first exceeded:",
    round(year_threshold_crossed, 2), "\n")


## 6.3 Maximum waiting time per calendar year ---------------------------------

waiting_time <- df_time_diag %>%
  dplyr::select(Start_blooddraw, dur_diag, value) %>%
  filter(!is.na(Start_blooddraw)) %>%
  mutate(calendar_year = round(Start_blooddraw / days_in_year + 2020, 0)) %>%
  group_by(calendar_year) %>%
  summarise(max_wait_days = max(dur_diag, na.rm = TRUE))

# Scale for dual-axis plot
y_primary_max   <- 6000
y_secondary_max <- 150
waiting_time <- waiting_time %>%
  mutate(
    scaled_wait = max_wait_days * (y_primary_max / y_secondary_max),
    QueueType   = "Max. waiting time"
  )


## 6.4 Outcomes: QALYs, life years, and costs ---------------------------------
# Number of analysis patients (burn-in years 1–5 and post-horizon years >15 excluded)
n_burnin   <- sum(arrival_times_trend <= 5 * days_in_year)
n_posthorizon <- sum(arrival_times_trend > 15 * days_in_year)
n_analysis <- length(arrival_times_trend) - n_burnin - n_posthorizon

# Exclude burn-in patients from all outcome calculations
df_attributes_filtered <- df_attributes %>%
  group_by(name) %>%
  filter(!any(time <= 5 * days_in_year)) %>%
  ungroup()

### QALYs
# average number of QALYs 
df_qaly <- df_attributes_filtered %>%
  filter(time > 5 * days_in_year, key == "QALY") %>%
  mutate(year = floor(time / days_in_year) - 5) %>%
  group_by(year) %>%
  summarise(total_q = sum(value), .groups = "drop")

qaly_per_patient <- sum(df_qaly$total_q) / n_analysis
qaly_per_patient

### Life years
# Average life years
df_ly <- df_attributes_filtered %>%
  filter(key %in% c("initialize_time", "PCa_death")) %>%
  dplyr::select(-value, -replication) %>%
  pivot_wider(names_from = key, values_from = time) %>%
  mutate(
    PCa_death  = pmin(coalesce(PCa_death, 15 * days_in_year), 15 * days_in_year),
    life_years = (PCa_death - initialize_time) / days_in_year
  ) %>%
  summarise(total = sum(life_years))

ly_per_patient <- df_ly$total / n_analysis
ly_per_patient

### Costs
# Average costs
df_cost <- df_attributes_filtered %>%
  filter(time > 5 * days_in_year & time <= 15 * days_in_year,
         grepl("Cost", key)) %>%
  mutate(year = floor(time / days_in_year) - 5) %>%
  group_by(year) %>%
  summarise(total_c = sum(value), .groups = "drop")

cost_per_patient <- sum(df_cost$total_c) / n_analysis
cost_per_patient


# 7. VISUALISATION ------------------------------------------------------------

okabe_ito_colors <- c(
  "Queue for Blood Draw"   = "#E69F00",
  "Queue for GP"           = "#56B4E9",
  "Queue for MRI"          = "#009E73",
  "Queue for Pathologist"  = "#F0E442",
  "Queue for Radiologist"  = "#0072B2",
  "Queue for Urologist"    = "#D55E00",
  "Max. waiting time"      = "#000000"
)

p_queue <- ggplot(df_queue, aes(x = (time / days_in_year) + 2020, y = value, color = key)) +
  geom_line(linewidth = 1.2) +
  geom_line(
    data        = waiting_time,
    aes(x = calendar_year, y = scaled_wait, color = QueueType),
    inherit.aes = FALSE,
    linewidth   = 1.2,
    linetype    = 1
  ) +
  scale_x_continuous(
    limits = c(2025, 2035),
    breaks = seq(2025, 2035, by = 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits   = c(0, y_primary_max),
    breaks   = seq(0, y_primary_max, by = 250),
    expand   = c(0, 0),
    sec.axis = sec_axis(~ . * (y_secondary_max / y_primary_max),
                        name = "Max. waiting time (days)")
  ) +
  scale_color_manual(values = okabe_ito_colors) +
  labs(
    title    = "Queue Size Over Time by Healthcare Resource",
    subtitle = "Simulation period 2025–2035",
    x        = "Year",
    y        = "Queue size (n)",
    color    = "Resource / metric"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(p_queue)