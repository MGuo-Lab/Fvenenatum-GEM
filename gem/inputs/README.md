# Inputs of the reconstruction

These files are the inputs of `curate_model.py` and `build_draft_baseline.py` that are not downloaded by
`download_data.py`. They are committed so that iJL1143 and the Draft model can be rebuilt without rerunning the
orthology step.

| File | Content | Produced by |
|---|---|---|
| `orthologs_cbs513_rbh_only.tsv` | Reciprocal best hits from blastp between the *A. niger* CBS 513.88 proteome and the *F. venenatum* A3/5 proteome (≥ 40 % identity, ≥ 70 % query coverage, E ≤ 1e-30): `template_gene`, `target_gene`, identity in both directions, E-value | `build_orthologs.py` (`jobs/orthology.sbatch`) |
| `Orthogroups.tsv` | OrthoFinder 2.5.5 orthogroups of the 10 proteomes listed in `download_data.py`, one column per species | OrthoFinder (`jobs/orthology.sbatch`) |
| `orthologs_cbs513.tsv` | The two combined, one row per *A. niger* gene with an *F. venenatum* ortholog; `method` is `RBH+OF` (both agree), `RBH`, `RBH(OF-disagree)` or `OF` (OrthoFinder pairwise orthologs only) | `merge_orthologs.py` |
| `fven_uniprot_ec.tsv` | UniProt entries of the A3/5 proteome UP000245910 with EC numbers, protein names and gene names, retrieved on 7 October 2026 from `https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000245910&fields=accession,ec,protein_name,gene_names&format=tsv` | UniProt |
| `secretion_whitelist.csv` | Metabolites allowed to leave the cell without a gene-supported transporter (Stage 3), with the basis and citation of each | curated by hand |
