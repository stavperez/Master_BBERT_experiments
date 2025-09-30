#!/bin/bash
#SBATCH --output=python_score_output_%A.txt  # Log file for output
#SBATCH --time=48:00:00                 # Set your desired time limit
#SBATCH --mem=216G                        # Memory requirements
#SBATCH --gres=gpu:1

DS=4
SIZE=128
# Process both R1 and R2
for READ_NUM in {1..2}; do
    # Construct input and output paths
#    INPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset2/output_reads_mix_256_orgs/output_reads_mix_256_orgs_R${READ_NUM}.fastq"
 #   OUTPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset2/bbert_results_256_orgs/bbert_scores_256_orgs_${READ_NUM}.csv"

#    INPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset_new/output_reads_mix_256_orgs/mix_${DS}_64bact_64euk_R${READ_NUM}.fastq"
 #   OUTPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset_new/output_reads_mix_256_orgs/bbert_results${DS}/bbert_scores_${DS}_${READ_NUM}.csv"

    INPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_${DS}_${SIZE}bact_${SIZE}euk_R${READ_NUM}.fastq"
    OUTPUT_FILE="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_size${SIZE}/bbert_scores_${SIZE}_${DS}_${READ_NUM}.csv"


    # Print the parameters for logging
    echo "Processing:"
    echo "Organism: $ORG"
    echo "Read Number: $READ_NUM"
    echo "Input File: $INPUT_FILE"
    echo "Output File: $OUTPUT_FILE"

    # Run the Python script
    srun python /sci/labs/aerez/stavperez/score.py "$INPUT_FILE" "$OUTPUT_FILE"

    # Check if the script executed successfully
    if [ $? -eq 0 ]; then
        echo "Processing completed successfully for ${ORG} R${READ_NUM}"
    else
        echo "Error occurred during processing for ${ORG} R${READ_NUM}"
        exit 1
    fi
done
