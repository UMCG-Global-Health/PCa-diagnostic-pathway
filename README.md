# Prostate Cancer Model

A discrete-event simulation (DES) model evaluating the **health and economic value of AI in the PCa diagnostic pathway** in the Dutch setting.

## Overview

This model simulates patient flow through a prostate cancer (PCa) diagnostic pathway, from GP referral through PSA/DRE testing, urologist assessment, MRI, biopsy, and treatment or active surveillance. It is implemented in R using the [`simmer`](https://r-simmer.org/) package.

The model compares four strategies, toggled via `mri_ai_reduction` and `rad_ai_reduction` in Section 1.6:

| Strategy | Description | `mri_ai_reduction` | `rad_ai_reduction` |
|----------|-------------|-------------------|-------------------|
| 0 | Base case (no AI) | 1 | 1 |
| 1 | AI-assisted MRI acquisition | 3 | 1 |
| 2 | AI-assisted radiologist reading | 1 | 2 |
| 3 | Combined AI-MRI + AI-reading | 3 | 2 |

AI effects are applied from simulation year 6 onward (years 1–5 serve as burn-in).

Key outcomes include discounted QALYs, life years, healthcare costs (all per analysis patient), resource queue lengths, and diagnostic waiting times.

---

## Model Structure

```
GP visit
  └─ PSA/DRE blood draw
       ├─ Negative → annual monitoring (up to 10 years)
       └─ Positive → Urologist intake
                       ├─ Negative → annual monitoring (up to 10 years)
                       └─ Positive → MRI scan + radiologist assessment
                                       ├─ Negative (PI-RADS 1–2) → monitor
                                       └─ Positive (PI-RADS ≥ 3) → Biopsy
                                                                       ├─ Negative → 2-year recheck
                                                                       └─ Positive → Treatment or Active Surveillance
```

Patient risk groups:

| Group | Description | Probability |
|-------|---------------------------------|-------------|
| 1 | Clinically significant PCa (Gleason ≥ 7) | 4.66% |
| 2 | Non-significant PCa | 6.44% |
| 3 | No prostate cancer | 88.9% |

---

## MRI Diagnostic Performance by Strategy

The MRI sensitivity and specificity differ across strategies. The values below reflect those reported in the manuscript (Table 2) and should be set in `mri_test_characteristics` (base case) and `mri_test_characteristics_ai` (AI strategies) in Section 1.3 of `model.R`.

### Strategy 0 — Base case (standard MRI + standard radiologist reading)

| Risk group | Sensitivity | Specificity |
|------------|-------------|-------------|
| Significant PCa | 0.91 | — |
| Non-significant PCa | 0.70 | — |
| No PCa | — | 0.58 |

- MRI scan time: 30 min (incl. 5 min preparation)
- Radiologist assessment time: 35 min

---

### Strategy 1 — AI-MRI (scan time reduction; −2% sensitivity and specificity vs. base case)

| Risk group | Sensitivity | Specificity |
|------------|-------------|-------------|
| Significant PCa | 0.89 | — |
| Non-significant PCa | 0.69 | — |
| No PCa | — | 0.57 |

- MRI scan time: 13 min (incl. 5 min preparation) → `mri_ai_reduction = 3`
- Radiologist assessment time: 35 min → `rad_ai_reduction = 1`

---

### Strategy 2 — AI-assessment (reading time reduction; +11.5 pp specificity, sensitivity unchanged vs. base case)

| Risk group | Sensitivity | Specificity |
|------------|-------------|-------------|
| Significant PCa | 0.91 | — |
| Non-significant PCa | 0.70 | — |
| No PCa | — | 0.69 |

- MRI scan time: 30 min → `mri_ai_reduction = 1`
- Radiologist assessment time: 17 min → `rad_ai_reduction = 2`

---

### Strategy 3 — Combined AI-MRI & AI-assessment (−2% sensitivity/specificity from AI-MRI + 11.5 pp specificity increase from AI-assessment)

| Risk group | Sensitivity | Specificity |
|------------|-------------|-------------|
| Significant PCa | 0.89 | — |
| Non-significant PCa | 0.69 | — |
| No PCa | — | 0.68 |

- MRI scan time: 13 min → `mri_ai_reduction = 3`
- Radiologist assessment time: 17 min → `rad_ai_reduction = 2`

> **Note:** The combined strategy specificity (0.68) reflects the net effect of the −2% reduction from AI-MRI applied on top of the +11.5 percentage point gain from AI-assessment relative to the base case (0.58 + 0.115 − 0.02 ≈ 0.68), as reported in Table 2 of the manuscript.

---

## Requirements

```r
install.packages(c("tidyverse", "simmer", "simmer.plot", "ggpubr", "here"))
```

Developed and tested with R 4.3+. Missing packages are installed automatically on first run.

---

## Usage

```r
source("R/model.R")
```

The script will:
1. Install any missing packages
2. Set all parameters and pre-sample patient risk groups
3. Build patient arrival schedules under the trend or flat volume scenario
4. Construct the simmer patient trajectory
5. Run the simulation (~16 years, 2020–2036)
6. Compute per-patient QALYs, life years, and costs for the analysis window (years 6–16)
7. Produce a queue-size plot with secondary axis for maximum diagnostic waiting time

---

### Switching strategies

To run a specific AI strategy, update three things in `model.R` before running:

1. **Section 1.3** — set `mri_test_characteristics_ai` to the sensitivity/specificity values for the chosen strategy (see table above)
2. **Section 1.6** — set `mri_ai_reduction` and `rad_ai_reduction` to the values for the chosen strategy
3. **Section 5** — choose the arrival scenario (`arrival_times_trend` or `arrival_times_flat`) in `add_generator()`

---

### Key parameters to adjust

| Parameter | Section | Description |
|-----------|---------|-------------|
| `mri_test_characteristics_ai` | 1.3 | MRI sensitivity/specificity under AI (strategy-dependent) |
| `mri_ai_reduction` | 1.6 | MRI scan time reduction factor (1 = no change, 3 = two-thirds reduction) |
| `rad_ai_reduction` | 1.6 | Radiologist reading time reduction factor (1 = no change, 2 = half) |
| `capacity_*` / `cap_*` | 1.6 | Resource capacities (GP, MRI, urologist, radiologist, pathologist) |
| `arrival_times_trend` / `arrival_times_flat` | 3 | Patient volume scenario |
| `discount_rate_qaly` / `discount_rate_cost` | 1.1 | Annual discount rates (%) |

---

## Outputs

| Object | Description |
|--------|-------------|
| `df_queue` | Queue length per resource over simulation time |
| `df_time_diag` | Per-patient diagnostic pathway duration |
| `waiting_time` | Maximum diagnostic waiting time per calendar year |
| `qaly_per_patient` | Mean discounted QALYs per analysis patient |
| `ly_per_patient` | Mean life years per analysis patient |
| `cost_per_patient` | Mean discounted costs (€) per analysis patient |
| `p_queue` | ggplot2 figure of queue sizes (2025–2035) |

### Analysis window

Patients entering in simulation years 1–5 (2020–2024) serve as burn-in to stabilise queues and are excluded from outcome calculations. Patients entering after year 15 are also excluded. The denominator `n_analysis` is derived directly from the arrival vector used in the simulation run.

---


## Citation

If you use this model in your work, please cite:

> van Dorst PWM, Fransen SJ, Vluttert T, Blanker MH, van Leeuwen PJ, Al-Uwini S, et al. The health and economic value of Artificial Intelligence in a resource-constrained diagnostic pathway of prostate cancer in the Netherlands – a discrete event simulation model.
---

## License

MIT License. See [LICENSE](LICENSE) for details.
