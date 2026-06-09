# Fig. 5 — Murine Validation

BXP-101, a multi-component NF-κB probe, modulates cognitive and molecular
phenotypes in an Aβ1–42 mouse model. (Used as a pharmacological probe, not a
PTGDS-specific intervention; no causal claim for PTGDS is made from these data.)

## Scripts

| File | Description |
|------|-------------|
| `Fig5_panels.R` | Generate individual panels A, B, C, E, F (no labels) |
| `Fig5_combine.R` | Assemble panels with labels into the final figure |

## Data

`data/murine/murine_behavioral_molecular_data.xlsx`

| Sheet | Panel | Contents |
|-------|-------|----------|
| Ab_Y-maze_summary | A | Spontaneous alternation (%) |
| Ab_PAT_summary | B | Passive avoidance retention latency |
| Ab_PCR_molecular | C, E, F | LCN2, TNF-α, IL-6, IL-1β, NGFR, BDNF, Keap1, APOE, ABCA1 |

> Panel D (LCN2/NGFR immunohistochemistry images) to be added upon completion.

## Usage
```r
Sys.setenv(MURINE_DATA = "data/murine/murine_behavioral_molecular_data.xlsx")
Sys.setenv(FIG5_OUT    = "output/Fig5")
source("analysis/Fig5_murine/Fig5_panels.R")
source("analysis/Fig5_murine/Fig5_combine.R")
```

## Groups

Sham / Aβ1–42 (20 μM, i.c.v.) / BXP-101 50/100/200/400 mg/kg / Donepezil 5 mg/kg
n = 8–10/group. One-way ANOVA + Tukey post-hoc.
