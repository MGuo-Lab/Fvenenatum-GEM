# Results

One folder per analysis. Each is written by the script of the same name in `analysis/scripts/` (the script name is
given where it differs) and holds a `provenance.json` with the commit of this repository, whether the working tree
had uncommitted changes, and the run parameters. The Draft model is labelled `draft` or `Draft`, and iJL1143
`curated`, `Curated` or `iJL1143`.

| Folder | Content | Used in the paper |
|---|---|---|
| `model_versions_qc` | Size, gene-association evidence, aerobic yields, maximum ATP yield and energy-generating-cycle tests of both models | Fig. 1, Fig. 2, Table 2 |
| `curation_stage_effects` | Each curation stage removed from iJL1143 or added to the Draft model: ATP yield, cycle tests, open exports, yields | Supplementary Table S3, Results |
| `energy_rules_effect` | Draft model with the directionality and carrier rules, run in the batch fermentations | Results |
| `flask_products` | Measured products, online data, yields and sugar depletion times per flask | Figs. 2, 4, Supplementary Fig. S2, Supplementary Table S5 |
| `flask_hplc` | Peak identity from the raw HPLC data: lactic acid window, 16.1 min peak areas, injection volumes | Fig. 4, Methods |
| `flask_dynamics` | dFBA of both models in the four fermentations in five settings; fitted initial biomass, final-biomass discrepancy, time-course RMSE and MAE, time courses | Fig. 3, Supplementary Table S5, Results |
| `lag_x0_fit` | Initial biomass of iJL1143 fitted per condition and shared (script `lag_x0_fit.py`) | Supplementary Fig. S1, Supplementary Note 1 |
| `methods_checks` | dFBA step size (0.5, 0.25, 0.1 h) and ash correction of the final-biomass comparison | Methods |
| `uptake_threshold_check` | Glucose threshold of sequential uptake (0.05, 0.5, 5 mM) and timing of infeasible dFBA steps | Methods |
| `overprediction` | Oxygen balance of the fermentations and P/O ratio sensitivity of iJL1143 (script `overprediction_tests.py`) | Methods, Results, Supplementary Note 2 |
| `yield_sensitivity` | Yield sensitivity to biomass composition and maintenance (±25 %) | Supplementary Fig. S2, Supplementary Note 2 |
| `km_sensitivity` | Sugar Km from 0.025 to 4.5 mM | Supplementary Note 2 |
| `glycerol_oxygen` | Glycerol timing against dissolved O2; O2-capped FBA and flux variability of both models (script `glycerol_oxygen_test.py`) | Results |
| `glycerol_diagnosis` | Net glycerol rates per sampling interval, their relation to dissolved O2, matched growth intervals, O2 shortfall test in iJL1143 | Results, Discussion |
| `peak_identity` | Retention and UV/RI ratio of the 16.1 min peak against single standards | Supplementary Fig. S3, Supplementary Note 3 |
| `ph_validation` | Medium pH from a charge balance driven by online biomass, four cumulative versions | Fig. 5, Table 1, Methods |
| `memote_ijl1143` | MEMOTE scores of both models (script `run_memote.py`) | Supplementary Table S4 |

`flask_hplc`, `glycerol_diagnosis` and `peak_identity` read the raw instrument files, which are not part of this
repository (see the main README).
