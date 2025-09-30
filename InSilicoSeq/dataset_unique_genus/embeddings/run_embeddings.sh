#!/bin/bash
#SBATCH --job-name=embeddings
#SBATCH --output=logs/embeddings_%j.txt
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH -t 7-0

set -euo pipefail

# ======= Directory config =======
BACT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact"
EUK_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Euk"
OUTPUT_BASE="output_scores"
PY_SCRIPT="BBERTooD/source/inference.py"

# ======= Read slice range =======
READ_START=0
READ_COUNT=2048

# ======= Organisms dictionary =======
declare -A organisms=(
#    [Actinomadura_madurae]="Bact"
 #   [Brevibacillus_brevis]="Bact"
  #  [Helicobacter_pylori]="Bact"
   # [Morganella_morganii]="Bact"
#    [Pseudomonas_aeruginosa]="Bact"
    [Escherichia_coli]="Bact"
    #[Acanthochromis_polyacanthus]="Euk"
    #[Saccharomyces_paradoxus]="Euk"
    #[Coffea_arabica]="Euk"
#    [Cyanidioschyzon_merolae_strain_10D]="Euk"
#    [Nakaseomyces_glabratus]="Euk"
    #[Ochotona_princeps]="Euk" # havent do yet
    #[Watersipora_subatra]="Euk"
)

# ======= Loop =======
for org_name in "${!organisms[@]}"; do
    type="${organisms[$org_name]}"
    
    if [[ "$type" == "Bact" ]]; then
        INPUT_DIR="$BACT_DIR"
    elif [[ "$type" == "Euk" ]]; then
        INPUT_DIR="$EUK_DIR"
    else
        echo "❌ Unknown type '$type' for $org_name - skipping."
        continue
    fi

    fq1="${INPUT_DIR}/${org_name}_R1.fastq"
    fq2="${INPUT_DIR}/${org_name}_R2.fastq"
    out_dir="${OUTPUT_BASE}/${org_name}"

    out1="${out_dir}/${org_name}_R1_scores.parquet"
    out2="${out_dir}/${org_name}_R2_scores.parquet"

    echo "- Processing: $org_name ($type)"
    echo "  Input files: $fq1 | $fq2"
    echo "  Output: $out1 | $out2"

    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "❌ Missing input files for $org_name - skipping."
        continue
    fi

    if [[ -f "$out1" && -f "$out2" ]]; then
        echo "✅ Already processed $org_name - skipping."
        continue
    fi

    mkdir -p "$out_dir"
    tmp_dir=$(mktemp -d)

    short_fq1="${tmp_dir}/${org_name}_R1_range.fastq"
    short_fq2="${tmp_dir}/${org_name}_R2_range.fastq"

    ## Compute line numbers (4 lines per read)
    start_line=$((READ_START * 4 + 1))
    end_line=$(( (READ_START + READ_COUNT) * 4 ))

    echo "  🔧 Taking reads $READ_START to $((READ_START + READ_COUNT - 1)) (lines $start_line to $end_line)"

    awk -v s=$start_line -v e=$end_line 'NR >= s && NR <= e' "$fq1" > "$short_fq1"
    awk -v s=$start_line -v e=$end_line 'NR >= s && NR <= e' "$fq2" > "$short_fq2"

    python "$PY_SCRIPT" \
        --input_dir "$tmp_dir" \
        --input_files "$(basename "$short_fq1")" "$(basename "$short_fq2")" \
        --output_dir "$out_dir" \
        --batch_size 2048 \
        --emb_out

    chmod -R 777 "$out_dir"
    rm -rf "$tmp_dir"

    echo "✅ Done with $org_name"
    echo
done
