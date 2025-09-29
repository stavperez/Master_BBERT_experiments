#!/bin/bash
#SBATCH --array=1-29
#SBATCH -o logs2/read_generator_%A_%a.out
#SBATCH -t 23:59:00
#SBATCH -c 4
#SBATCH --mem=23999

echo "Starting job with SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

#LINE_NUM=$((SLURM_ARRAY_TASK_ID))
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
#DATASET_NUM=4

#CSV_FILE="ds${DATASET_NUM}_shuffled_bact_org_info.csv"  # or .bact.
#CSV_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset_embeddings/ds_embeddings.csv"

#CSV_FILE="/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv"
#CSV_FILE="/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus.csv"
CSV_FILE="/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus_missing_or_invalid.csv"

#CSV_FILE="chromosome_assemblies_unique_genus_with_paths.csv"
#CSV_FILE="/cs/usr/stavperez/sp/genomes_new/genome_refseq_eukaryotes.csv"
#TEMP_CSV="/cs/usr/stavperez/sp/InSilicoSeq/dataset${DATASET_NUM}/tmp/temp_row_${SLURM_ARRAY_TASK_ID}.csv"
TEMP_CSV="/cs/usr/stavperez/sp/genomes_new/tmp/temp1_row_${SLURM_ARRAY_TASK_ID}.csv"
OUTPUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/"

# Extract header + row
head -n 1 "$CSV_FILE" > "$TEMP_CSV"
sed -n "${LINE_NUM}p" "$CSV_FILE" >> "$TEMP_CSV"

if [ ! -s "$TEMP_CSV" ]; then
    echo "Error: No data found for row ${LINE_NUM}. Exiting."
    exit 1
fi

# Infer organism type from file name
if [[ "$CSV_FILE" == *"euk"* ]]; then
    ORG_TYPE="Euk"
elif [[ "$CSV_FILE" == *"bact"* ]]; then
    ORG_TYPE="Bact"
else
#    echo "Could not infer organism type from file name."
#    exit 1
    ORG_TYPE="Bact"
fi

# Run the Python script
cmd="python -u read_generator_one_org.py --csv \"$TEMP_CSV\" --org \"$ORG_TYPE\" --output_dir \"$OUTPUT_DIR\""
echo "Executing: $cmd"
eval $cmd

# Clean up
rm -f "$TEMP_CSV"
echo "Job completed successfully for row ${LINE_NUM}."
