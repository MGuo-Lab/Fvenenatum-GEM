#!/bin/bash
#SBATCH --job-name=fven_orthology
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --output=logs/orthology-%j.out
#SBATCH --error=logs/orthology-%j.err

# Orthology between the template (CBS 513.88 version of iJB1325) and F. venenatum A3/5:
# reciprocal best hits from blastp (40 % identity, 70 % query coverage, E <= 1e-30) and OrthoFinder orthogroups
# across 10 fungal proteomes, merged by merge_orthologs.py. Run from gem/ after download_data.py:
#   sbatch jobs/orthology.sbatch
# Environment used for the published tables: BLAST+ 2.16.0, OrthoFinder 2.5.5, DIAMOND 2.2.8, Python 3.11 with
# COBRApy 0.32.1, numpy 1.26.4 and scipy 1.11.3 (OrthoFinder 2.5.5 fails on numpy 2 in its orthologue step).
# The resulting tables are committed in inputs/ (orthologs_cbs513.tsv, orthologs_cbs513_rbh_only.tsv and the
# Orthogroups.tsv of OrthoFinder), so the later steps do not need this job.

set -euo pipefail
mkdir -p logs orthology
T=${SLURM_CPUS_PER_TASK:-16}

echo "[$(date)] RBH blastp on $T cores"
python build_orthologs.py --template-model raw_data/MOESM3_cbs51388.xml \
    --aniger-fasta raw_data/aniger_cbs513.faa \
    --output orthology/orthologs_cbs513_rbh.tsv --threads "$T" --keep-blast

echo "[$(date)] OrthoFinder"
rm -rf orthology/orthofinder_output
orthofinder -f raw_data/orthofinder_input -t "$T" -a 8 -o orthology/orthofinder_output

echo "[$(date)] merge"
python merge_orthologs.py --rbh-tsv orthology/orthologs_cbs513_rbh.tsv --of-dir orthology/orthofinder_output \
    --aniger-species Aniger_CBS513 --output orthology/orthologs_cbs513.tsv \
    --backup orthology/orthologs_cbs513_rbh_only.tsv
echo "[$(date)] done"
