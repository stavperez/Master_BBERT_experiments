#!/bin/bash
#SBATCH --job-name=spades_assembly
#SBATCH --output=spased_%j    #spades_log_real_data_%j.out
#SBATCH --time=98:00:00       # Job runtime
#SBATCH --mem=640G           #128G              # Memory
#SBATCH --cpus-per-task=80    #8      # Number of CPUs

# Load SPAdes environment
#module load python3

DS=4
SIZE=128
  # Define input and output file paths based on the organism
PREFIX="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/"
#R1_FILE="${PREFIX}/bbert_results${DS}_size${SIZE}/classification_output_bacteria_R1.fastq"
#R2_FILE="${PREFIX}/bbert_results${DS}_size${SIZE}/classification_output_bacteria_R2.fastq"

R1_FILE="${PREFIX}/bbert_results${DS}_size${SIZE}/new_classification_output_bacteria_R1.fastq"
R2_FILE="${PREFIX}/bbert_results${DS}_size${SIZE}/new_classification_output_bacteria_R2.fastq"

#R1_FILE="${PREFIX}/mix_${DS}_${SIZE}bact_${SIZE}euk_R1.fastq"
#R2_FILE="${PREFIX}/mix_${DS}_${SIZE}bact_${SIZE}euk_R2.fastq"

#OUTPUT_DIR="${PREFIX}/assembly_output/classified_bact_${DS}_${SIZE}bact_${SIZE}euk"
OUTPUT_DIR="${PREFIX}/assembly_output/new_classified_bact_${DS}_${SIZE}bact_${SIZE}euk"
#OUTPUT_DIR="${PREFIX}/assembly_output/mix_${DS}_${SIZE}bact_${SIZE}euk"


  # Run SPAdes
python /cs/usr/stavperez/sp/SPAdes-4.0.0-Linux/bin/spades.py \
  -1 "$R1_FILE" \
  -2 "$R2_FILE" \
  -o "$OUTPUT_DIR" \
  -t 80 \
  -m 640 \
  --meta

