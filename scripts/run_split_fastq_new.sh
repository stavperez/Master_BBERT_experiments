#!/bin/bash
DS="4"
SIZE="128"
R1_FASTQ="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_${DS}_${SIZE}bact_${SIZE}euk_R1.fastq"
R2_FASTQ="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_${DS}_${SIZE}bact_${SIZE}euk_R2.fastq"
SCORE_CSV_R1="/sci/backup/aerez/aerez/from_Dmitry/scores/bact_frame_scores/mix_${DS}_${SIZE}bact_${SIZE}euk_R1_scores.csv"
SCORE_CSV_R2="/sci/backup/aerez/aerez/from_Dmitry/scores/bact_frame_scores/mix_${DS}_${SIZE}bact_${SIZE}euk_R2_scores.csv"
OUTPUT_PREFIX="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_size${SIZE}/new_classification_output"

# Run the Python script
python3 split_fastq_by_score_new.py \
    --r1 "$R1_FASTQ" \
    --r2 "$R2_FASTQ" \
    --csv_r1 "$SCORE_CSV_R1" \
    --csv_r2 "$SCORE_CSV_R2" \
    --output "$OUTPUT_PREFIX"

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "FASTQ splitting completed successfully!"
else
    echo "Error: FASTQ splitting failed!" >&2
    exit 1
fi
