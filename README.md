# Fvenenatum-GEM

iJL1143, a genome-scale metabolic model of the mycoprotein fungus *Fusarium venenatum* A3/5, together with the
code, data and results of the accompanying paper:

> Lee J., Vinestock T., Guo M. Rule-based curation and validation in batch fermentation of a genome-scale metabolic
> model of the mycoprotein fungus *Fusarium venenatum* A3/5.

iJL1143 was reconstructed from the CBS 513.88 version of the *Aspergillus niger* model iJB1325
(Brandl et al. 2018) by orthology and curated with six stages of rules based on homology, thermodynamics and the
literature. The paper compares it with the uncurated Draft model in four oxygen-enriched batch fermentations of
strain A3/5 on glucose and fructose (Vinestock & Guo 2026).

## Contents

| Folder | Contents |
|---|---|
| `model/` | `iJL1143.xml` (Curated model) and `Draft.xml` (Draft model), SBML Level 3 with FBC version 2 |
| `model/curation/` | Curation change log and the checks behind the gene layer of iJL1143 |
| `data/flask/` | Batch fermentation data: HPLC time courses, online biomass, pH and dissolved O2, end-point dry weight |
| `results/` | Result files behind every figure, table and number in the paper, one folder per analysis |
| `figures/` | Figures of the paper (PDF, SVG, PNG) and CSV versions of Tables 1, 2 and S1 |
| `gem/` | Reconstruction, curation and dFBA code (Python package `fven-gem`) and its committed inputs |
| `analysis/` | Data loading, simulations, analyses and figures (Python package `fbaval`) |

## The models

| | Draft | Curated (iJL1143) |
|---|---|---|
| Reactions | 2,321 | 2,264 |
| Metabolites | 1,818 | 1,797 |
| Genes | 1,034 | 1,143 |
| Exchange reactions open for secretion | 377 | 47 |
| Maximum ATP yield (mol ATP per mol glucose) | 49.1 | 29.9 |
| Aerobic biomass yield on glucose (g g⁻¹) | 0.599 | 0.522 |

- **iJL1143** (`model/iJL1143.xml`, version 1.1.0) contains all six curation stages. Its genes are UniProt
  accessions of the A3/5 proteome (UniProt UP000245910; genome ENA PRJEB7533, King et al. 2018), written as
  `G_<accession>` in the SBML file.
- **Draft** (`model/Draft.xml`) keeps every reaction of the template, with gene associations rewritten through
  reciprocal best hits only. It has the biomass composition and maintenance of iJL1143, so the two models differ only
  in their network and gene associations.
- Both models use the reaction and metabolite identifiers of iJB1325. Exchange reactions are named
  `BOUNDARY_<metabolite>e`, for example `BOUNDARY_GLCe` (glucose), `BOUNDARY_FRUe` (fructose) and `BOUNDARY_O2e`
  (oxygen). Biomass is formed by `GROWTH` and leaves through `DRAIN_Biomass`, the objective. Biomass is normalised
  to 1 g mmol⁻¹, so the objective value is the specific growth rate in h⁻¹. Growth-associated maintenance is
  61 mmol ATP gDW⁻¹ in `GROWTH`; non-growth-associated maintenance is the reaction `NGAM`, fixed at
  1.9 mmol ATP gDW⁻¹ h⁻¹.
- The default medium is aerobic minimal glucose medium: glucose uptake up to 10 mmol gDW⁻¹ h⁻¹, with ammonium,
  phosphate, sulfate and O2.

```python
import cobra

model = cobra.io.read_sbml_model("model/iJL1143.xml")
model.reactions.BOUNDARY_GLCe.lower_bound = -1.0      # glucose uptake, mmol gDW-1 h-1
mu = model.slim_optimize()                             # 0.094 h-1
print(mu / 0.18016)                                    # biomass yield, 0.522 g per g glucose
```

### Curation change log

`model/curation/curation_log.csv` lists every change from the template to iJL1143, one row per reaction and rule,
with the columns `reaction_id`, `stage`, `action`, `rule`, `evidence` and `confidence` (Thiele & Palsson 2010
scores: 4 biochemical, 3 genetic, 2 sequence or physiological, 1 modelling only, 0 none).

| `stage` | Stage in the paper |
|---|---|
| `1_draft` | 1, orthology filter (tier 1 reciprocal best hit, tier 2 OrthoFinder orthogroup) |
| `2_ec_rescue`, `2_gapfill` | 2, reactions restored by UniProt EC annotation and gap-filling for biomass synthesis |
| `3_secretion` | 3, secretion only through a gene-supported transporter or with literature support |
| `4_redox`, `4_thermo` | 4, redox rules and directionality of NAD(P)-dependent reactions |
| `5_biomass` | 5, biomass composition of growing mycelium and NGAM |
| `6_v2_thermo`, `6_v2_acid_export`, `6_v3_gtp_carrier` | 6, directionality of all metabolic reactions, ATP cost of organic-acid export, direction of the GTP/GDP carrier |

The other files in `model/curation/`:

- `draft_evidence.tsv`: orthology evidence for each template reaction in Stage 1.
- `energy_leaks.tsv`: the cofactor-route tests of Stage 4.
- `unmapped_subunit_hits.csv`: Smith–Waterman alignments of the nine complex subunits without an ortholog.
- `release_gene_changes.csv`: the 61 gene objects that the release step removed or replaced. Only the gene layer
  changes in this step; reactions, bounds and metabolites stay as they are.

## Batch fermentation data

`data/flask/` holds the four batch fermentations of Vinestock & Guo (2026), Glc 15, Glc 10 + Fru 5, Glc 5 + Fru 10
and Fru 15 (g L⁻¹), with three flasks each. [`data/README.md`](data/README.md) describes the culture conditions, the
files and their units.

## Results, figures and tables

Each folder in `results/` is written by the script of the same name in `analysis/scripts/` and holds a
`provenance.json` with the commit of this repository that produced it. [`results/README.md`](results/README.md)
lists the folder behind each figure, table and statement of the paper.

| Paper | File in `figures/` | Results | Script |
|---|---|---|---|
| Fig. 1 | `fig1_reconstruction` | `model_versions_qc`, model files, change log | `paper_figures.py` |
| Fig. 2 | `fig2_curation_effects` | `model_versions_qc`, `flask_products` | `paper_figures.py` |
| Fig. 3 | `fig3_fermentation_time_courses` | `flask_dynamics` | `paper_figures.py` |
| Fig. 4 | `fig4_products` | `flask_products`, `flask_hplc` | `paper_figures.py` |
| Fig. 5 | `fig5_ph` | `ph_validation` | `paper_figures.py` |
| Supplementary Fig. S1 | `figS1_lag` | `lag_x0_fit` | `lag_x0_figure.py` |
| Supplementary Fig. S2 | `figS2_yield_sensitivity` | `yield_sensitivity`, `flask_products` | `paper_figures.py` |
| Supplementary Fig. S3 | `figS3_arabitol_identity` | `peak_identity` | `paper_figures.py` |
| Table 1 | `tables/table1_inputs.csv` | `ph_validation` (medium and ionic strength) | `paper_figures.py` |
| Table 2 | `tables/table2_models.csv` | `model_versions_qc`, model files | `paper_figures.py` |
| Supplementary Table S1 | `tables/tableS1_curation_rules.csv` | change log | `paper_figures.py` |
| Supplementary Table S2 | | `BIOMASS_COMPOSITION` and `COMPOSITION_SOURCE` in `gem/curate_model.py` | |
| Supplementary Table S3 | | `curation_stage_effects/effects.csv` | `curation_stage_effects.py` |
| Supplementary Table S4 | | `memote_ijl1143/scores.csv` | `run_memote.py` |
| Supplementary Table S5 | | `flask_dynamics/metrics.csv`, `flask_products/replicates.csv` | `flask_dynamics.py`, `flask_products.py` |

The files in `figures/` are the figures of the manuscript.

## Reproducing the results

The results were computed with Python 3.12 and the package versions in [`requirements.txt`](requirements.txt).
With [uv](https://docs.astral.sh/uv/):

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
uv pip install --python .venv/bin/python -e gem -e analysis
source .venv/bin/activate
```

`python -m venv .venv` followed by `pip install` with the same arguments works as well. COBRApy installs GLPK
(`swiglpk`).

**Template and proteomes.** `python gem/download_data.py` downloads the template model and the 10 proteomes into
`gem/raw_data/` and checks each against the version used for iJL1143. `--only template` fetches the template alone,
which is enough for everything except the orthology step.

**Models.** The orthology tables are committed in `gem/inputs/`, so the models can be rebuilt without BLAST or
OrthoFinder:

```bash
cd gem
python curate_model.py              # stages 1-6 -> build/curated.xml and model/curation/
python check_unmapped_subunits.py   # Smith-Waterman check of the nine unmapped subunits (Biopython)
python build_release_model.py       # gene layer of the release -> model/iJL1143.xml
python build_draft_baseline.py      # -> model/Draft.xml
```

The rebuilt models agree with the committed ones in every reaction, bound, stoichiometric coefficient, gene
association, metabolite and gene. Curation computes the reversibility index with the eQuilibrator API, which
downloads its compound cache (about 1.3 GB) on first use. The orthology step (`gem/jobs/orthology.sbatch`:
reciprocal best hits from blastp, OrthoFinder across 10 proteomes and `merge_orthologs.py`) needs BLAST+,
OrthoFinder and DIAMOND and was run on a computing cluster.

**Analyses and figures.** `bash analysis/run_analyses.sh` runs every analysis in dependency order and then draws the
figures; on 60 cores it takes about 20 min. Single analyses run as `python scripts/<name>.py` from `analysis/`.
Three analyses (`flask_hplc`, `glycerol_diagnosis` and `peak_identity`) read the raw instrument files, that is the
Shimadzu LabSolutions HPLC data files and the PreSens exports. These files (about 1 GB) are not part of the
repository; `run_analyses.sh` skips the three steps unless `FBAVAL_RAW_DIR` points to them, and the later steps then
use the committed result files of these steps. Tests run with `pytest` in `gem/` and in `analysis/`.

## Software

| Software | Version | Role |
|---|---|---|
| BLAST+ | 2.16.0 | Reciprocal best hits (blastp) |
| OrthoFinder, DIAMOND | 2.5.5, 2.2.8 | Orthogroups across 10 fungal proteomes |
| Biopython | 1.88 | Smith–Waterman alignment of unmapped complex subunits |
| eQuilibrator API | 0.8.1 | Reversibility index for the directionality rules |
| COBRApy, GLPK | 0.32.1, 5.0 | Model input and output, flux balance analysis, dFBA |
| MEMOTE | 0.17.0 | Model quality scores |
| SciPy, NumPy, pandas, Matplotlib | 1.18.1, 2.5.3, 2.3.3, 3.11.2 | Fitting, data handling, figures |
| Python | 3.12.3 (orthology step 3.11) | |

## References

- Brandl J., Aguilar-Pontes M.V., Schäpe P., Noerregaard A., Arvas M., Ram A.F.J., Meyer V., Tsang A., de Vries
  R.P., Andersen M.R. (2018) A community-driven reconstruction of the *Aspergillus niger* metabolic network.
  *Fungal Biol. Biotechnol.* 5, 16. https://doi.org/10.1186/s40694-018-0060-7
- King R., Brown N.A., Urban M., Hammond-Kosack K.E. (2018) Inter-genome comparison of the Quorn fungus *Fusarium
  venenatum* and the closely related plant infecting pathogen *Fusarium graminearum*. *BMC Genomics* 19, 269.
  https://doi.org/10.1186/s12864-018-4612-2
- Thiele I., Palsson B.Ø. (2010) A protocol for generating a high-quality genome-scale metabolic reconstruction.
  *Nat. Protoc.* 5, 93–121. https://doi.org/10.1038/nprot.2009.203
- Vinestock T., Guo M. (2026) Experiments & modelling of batch fermentation of *Fusarium venenatum* on
  glucose-fructose mixtures. *Syst. Control Trans.* 5, 776–782. https://doi.org/10.69997/sct.178293

## Contact

Jaewook Lee and Miao Guo, Department of Engineering, King's College London.
