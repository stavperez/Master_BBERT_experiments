#!/bin/bash
#SBATCH --job-name=embeddings_confusion
#SBATCH --output=logs/embeddings_confusion_%j.txt
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH -t 7-0

set -euo pipefail

# ======= Input files =======
#INPUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings"
#FQ1="reads_by_confusion_type_R1.fastq"
#FQ2="reads_by_confusion_type_R2.fastq"
INPUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/reads_by_reading_frame"
FQ1="reading_frame_-2.0_R1.fastq"
FQ2="reading_frame_-2_R2.fastq"

# ======= Output config =======
OUTPUT_DIR="${INPUT_DIR}/output_scores_embeddings"
PY_SCRIPT="/sci/labs/aerez/alekhin_dm_81/projects/BBERTooD/source/inference.py"

# ======= Create output directory =======
mkdir -p "$OUTPUT_DIR"

echo "- Running embeddings on confusion-type FASTQs"
echo "  Input files: $INPUT_DIR/$FQ1 and $FQ2"
echo "  Output dir: $OUTPUT_DIR"

# ======= Run inference =======
python "$PY_SCRIPT" \
    --input_dir "$INPUT_DIR" \
    --input_files "$FQ1" "$FQ2" \
    --output_dir "$OUTPUT_DIR" \
    --batch_size 2048 \
    --emb_out

chmod -R 777 "$OUTPUT_DIR"

echo "✅ Finished embedding generation"
