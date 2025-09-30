#!/bin/bash
#SBATCH --job-name=align_contigs_adj_errors
#SBATCH --output=align_contigs_adj_errors_%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=32G

DS=4
SIZE=128
BLAST_DB="/cs/usr/stavperez/sp/genomes_new/blastdb/genome_unique_genus_bacts_master_db"

#BLAST_DB="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_${SIZE}/genome_>
THREADS=8
EVALUE_THRESHOLD="1e-5"

declare -A CONTIG_SETS=(
  ["new_classified_bact_${DS}_${SIZE}bact_${SIZE}euk"]="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/assembly_output/new_classified_bact_${DS}_${SIZE}bact_${SIZE}euk/contigs.fasta"
  #["mix_${DS}_${SIZE}bact_${SIZE}euk"]="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/assembly_output/mix_${DS}_${SIZE}bact_${SIZE}euk/contigs.fasta"
)

OUT_DIR="/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results${DS}_${SIZE}/blast_results_adj_errors"
mkdir -p "$OUT_DIR"


for name in "${!CONTIG_SETS[@]}"; do
  fasta="${CONTIG_SETS[$name]}"
  raw_out="${OUT_DIR}/${name}_blast_raw.tsv"
  final_csv="${OUT_DIR}/${name}_blast_with_error_types.csv"

  echo "Running BLAST for: $name"
  blastn \
    -query "$fasta" \
    -db "$BLAST_DB" \
    -out "$raw_out" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -max_target_seqs 10 \
    -evalue $EVALUE_THRESHOLD \
    -num_threads $THREADS

  echo "Post-processing error types in: $raw_out"
  python3 - <<EOF
import pandas as pd

def classify_alignment_errors(qseq, sseq):
    single_errors = 0
    adjacent_pairs = 0
    longer_runs = 0
    run = 0

    for q, s in zip(qseq, sseq):
        if q != s or q == '-' or s == '-':
            run += 1
        else:
            if run == 1:
                single_errors += 1
            elif run == 2:
                adjacent_pairs += 1
            elif run > 2:
                longer_runs += 1
            run = 0

    if run == 1:
        single_errors += 1
    elif run == 2:
        adjacent_pairs += 1
    elif run > 2:
        longer_runs += 1

    return pd.Series([single_errors, adjacent_pairs, longer_runs])

df = pd.read_csv("$raw_out", sep="\t", header=None)
df.columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq"
]

df[["Single_Errors", "Adjacent_Pairs", "Longer_Runs"]] = df.apply(
    lambda row: classify_alignment_errors(row["qseq"], row["sseq"]), axis=1)

df.drop(columns=["qseq", "sseq"], inplace=True)
df.to_csv("$final_csv", index=False)
print(f"Saved: $final_csv")
EOF

done

echo "All done."
