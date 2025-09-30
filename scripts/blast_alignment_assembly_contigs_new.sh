#!/bin/bash
#SBATCH --job-name=align_contigs
#SBATCH --output=align_contigs_%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=32G


#### Align top 10 contigs ####

# === CONFIG ===
DS=4
SIZE=128
BLAST_DB="/cs/usr/stavperez/sp/genomes_new/blastdb/genome_unique_genus_bacts_master_db"

#BLAST_DB="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_${SIZE}/genome_bacts_db_${DS}_${SIZE}"
THREADS=8
#EVALUE_THRESHOLD="1e-10"
EVALUE_THRESHOLD="1e-3"
OUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_${SIZE}/blast_results_bacts_db"
mkdir -p "$OUT_DIR"

declare -A CONTIG_SETS=(
  ["new_classified_bact_${DS}_${SIZE}bact_${SIZE}euk"]="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/assembly_output/new_classified_bact_${DS}_${SIZE}bact_${SIZE}euk/contigs.fasta"
#  ["mix_${DS}_${SIZE}bact_${SIZE}euk"]="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/assembly_output/mix_${DS}_${SIZE}bact_${SIZE}euk/contigs.fasta"
)

for name in "${!CONTIG_SETS[@]}"; do
  fasta="${CONTIG_SETS[$name]}"
  raw_out="${OUT_DIR}/${name}_blast_raw.tsv"
  filtered_csv="${OUT_DIR}/${name}_blast_final.csv"

  echo "Running BLAST for: $name"
  blastn \
    -query "$fasta" \
    -db "$BLAST_DB" \
    -out "$raw_out" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 10 \
    -evalue 1e-5 \
    -num_threads "$THREADS"

  echo "Filtering by e-value < $EVALUE_THRESHOLD and saving as CSV..."

  python3 - <<EOF
import pandas as pd

tsv_file = "$raw_out"
output_csv = "$filtered_csv"
evalue_thresh = float("$EVALUE_THRESHOLD")

columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
]

df = pd.read_csv(tsv_file, sep="\t", header=None, names=columns)
df = df[df["evalue"] < evalue_thresh]
df.to_csv(output_csv, index=False)
print(f"Saved filtered hits to: {output_csv}")
EOF

done

echo "All done."
