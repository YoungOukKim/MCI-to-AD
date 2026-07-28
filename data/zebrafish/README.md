# Zebrafish Experimental Data

Experimental data from the zebrafish MCI model (Zefit Inc., Daegu, Republic of Korea).
All experiments performed at 14 dpf and 21 dpf using *Danio rerio* larvae.

---

## MCI Induction Protocol

| Agent | Concentration | Exposure window |
|---|---|---|
| LPS (*E. coli* O111:B4, Sigma-Aldrich Cat# L2630) | 5 μg/L | 2–14 dpf |
| D-galactose (Sigma-Aldrich Cat# G0750) | 0.2 mg/L | 4–14 dpf |
| High-fat diet (HFD; Gemma micro ZF 75 + 40% egg yolk) | 10% | 4–14 dpf |

**BXP-101 composition** (co-treatment 6–14 dpf):

| Dose | Atractylodin | Wedelolactone | Honokiol | Total |
|---|---|---|---|---|
| 0.3 μg/ml | 0.125 μg/ml | 0.125 μg/ml | 0.050 μg/ml | 0.30 μg/ml |
| 0.4 μg/ml | 0.167 μg/ml | 0.167 μg/ml | 0.067 μg/ml | 0.40 μg/ml |
| 0.6 μg/ml | 0.250 μg/ml | 0.250 μg/ml | 0.100 μg/ml | 0.60 μg/ml |

---

## File Descriptions

### behavior_14dpf.csv
- **Assay**: Red ball avoidance (visual stimuli avoidance response)
- **Method**: 6-well plate, 5 larvae/well, 4 wells per biological replicate = 20 larvae per group per replicate; 3 biological replicates
- **Measurement**: % of larvae in non-stimuli area (6 min, 12-sec intervals, mean)
- **Groups**: Control / MCI / MCI+Donepezil(10 μM) / MCI+BXP101(0.3/0.4/0.6 μg/ml)
- **Replicates**: 3 biological replicates × 4 wells each
- **Statistics**: Two-tailed t-test; `*` p<0.05 vs Control; `#` p<0.05 vs MCI

### qPCR_14dpf.csv
- **Assay**: RT-qPCR gene expression at 14 dpf
- **Method**: ddCt method, relative expression vs Control (set to 1.000)
- **Genes**: ptgdsb.1, ptgdsb.2, ngfr, bdnf, gfap, tnf_alpha, il1b, il6, cox2_ptgs2a, creb, dlg4_psd95 (psd-95), nestin
- **Groups**: Control / MCI / MCI+Donepezil(10 μM) / MCI+BXP101(0.3/0.4/0.6 μg/ml)
- **Replicates**: 3 biological replicates per group
- **Statistics**: Two-tailed t-test; `*` p<0.05, `**` p<0.01 vs Control; `#` p<0.05, `##` p<0.01 vs MCI

### qPCR_21dpf.csv
- **Assay**: RT-qPCR gene expression at 21 dpf
- **Method**: ddCt method, relative expression vs Control
- **Groups**: Control / MCI / MCI+BXP101(0.4 μg/ml) only
- **Replicates**: 4 biological replicates per group
- **Key finding**: ptgdsb.1/2 in MCI rises across the larval window — modestly below Control at 14 dpf (0.87/0.85) and slightly above Control by 21 dpf (1.063/1.066) — consistent with an ascending compensatory limb rather than a steady-state elevation. See `qPCR_14dpf.csv` for the 14 dpf values.
- **BXP-101 lowers larval ptgdsb at 14 dpf in a graded manner** (ptgdsb.1: 0.96 / 0.72 / 0.65 at 0.3 / 0.4 / 0.6 μg/ml), the direction expected if reducing inflammatory load reduces compensatory demand on PTGDS rather than if the compound acted on PTGDS directly

### ELISA_14dpf.csv
- **Assay**: Sandwich ELISA for cytokine protein quantification at 14 dpf
- **Kits**: Zebrafish TNFα ELISA Kit (ELK8512-96T) / Zebrafish IL-6 ELISA Kit (ELK9663-96T), ELK Biotechnology
- **Sample**: 30 larvae pooled per group; PBS homogenate w:v = 1:9. **Note:** a single pooled homogenate per group, i.e. no biological replicates — ELISA values are descriptive magnitudes and no statistics are computed from them
- **Unit**: pg/ml
- **Key finding**: TNF-α protein is suppressed by BXP-101 at 0.4 and 0.6 μg/ml (0.67× and 0.55× MCI) while TNF-α mRNA at 0.6 μg/ml is unchanged relative to MCI (1.28 vs 1.28 vs Control; ratio 1.00× MCI), indicating post-transcriptional rather than transcriptional suppression. **Correction (2026-07):** an earlier version of this note stated "TNF-α mRNA (BXP0.6: 1.28× MCI)"; 1.28 is the value relative to Control, not to MCI. The mRNA ratio to MCI is 1.00×.
- **Dose response is threshold-type, not monotonic**: at 0.3 μg/ml TNF-α (41.9 pg/ml) and IL-6 (94.7 pg/ml) were marginally above MCI (39.7 and 91.8); suppression appears only at 0.4 and 0.6 μg/ml

### WIF_14dpf.csv
- **Assay**: Whole-mount double immunofluorescence at 14 dpf
- **Markers**: BLBP (radial glia / reactive gliosis marker) + Nestin (neural progenitor marker)
- **Antibodies**: anti-BLBP (Abcam ab32423, TRITC secondary ab6718); anti-Nestin clone 10C2 (Merck MAB5326, Alexa Fluor 488 ab150113)
- **Imaging**: Agilent BioTek Lionheart FX automated microscope, identical settings across groups
- **Quantification**: ImageJ, fixed ROI (whole head region), mean gray value (arbitrary units, a.u.)
- **N**: 10 larvae per group
- **Statistics**: Two-tailed t-test; `**` p<0.01, `***` p<0.001 vs Control; `###` p<0.001 vs MCI

---

## Revision note (2026-07)

This README was corrected during manuscript revision. Changes: (1) the TNF-α mRNA ratio in the ELISA section was misstated (1.28 is relative to Control, not to MCI); (2) the 21 dpf key-finding now states the 14 dpf values as well, so the ascending trajectory is visible from the README alone; (3) the threshold-type (non-monotonic) dose response at 0.3 μg/ml is stated explicitly; (4) the ELISA is identified as a single pooled measurement per group without replicates; (5) behavioural sample size is stated unambiguously. No data files were altered.
