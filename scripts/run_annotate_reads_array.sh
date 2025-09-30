#!/bin/bash
#SBATCH --array=700-800
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=logs1/annotate_reads_%A_%a.out


# File paths
READS_CSV="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_1_128bact_128euk_R1.csv"
BACT_CSV="/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv"
EUK_CSV="/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus.csv"
#OUTPUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/annotated_reads_split_mix_1"
OUTPUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/annotated_reads_split_mix_1_non_coding"
#SCRIPT="annotate_single_organism.py"
SCRIPT="annotate_non_coding_single_organism.py"

mkdir -p $OUTPUT_DIR

# Get organism name for this SLURM task
ORG_LIST=$(python -c "
import pandas as pd
b = pd.read_csv('$BACT_CSV')
e = pd.read_csv('$EUK_CSV')
df = pd.concat([b, e])
print(df['organism_name'].iloc[$SLURM_ARRAY_TASK_ID - 1])
")
echo "Processing organism: $ORG_LIST"

# Run annotation for this organism
python $SCRIPT --reads $READS_CSV \
               --organism "$ORG_LIST" \
               --bacteria_csv $BACT_CSV \
               --eukaryote_csv $EUK_CSV \
               --output_dir $OUTPUT_DIR \
               --skip_existing

