#!/bin/bash

# Exit if anything fails
set -e

# === CONFIG ===
BACT_CSV="/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv"
COMBINED_FNA="/cs/usr/stavperez/sp/genomes_new/combined_unique_genus_bact_genomes.fna"
BLAST_DB="/cs/usr/stavperez/sp/genomes_new/blastdb/genome_unique_genus_bacts_master_db"
mkdir -p "$(dirname "$BLAST_DB")"

# === STEP 1: Extract .fna.gz paths for the first 256 bacteria ===
echo "Combining .fna.gz files from first 256 bacteria in CSV..."

# Remove previous combined file
rm -f "$COMBINED_FNA"

# Extract paths and concatenate genomes
python3 <<EOF
import pandas as pd
df = pd.read_csv("$BACT_CSV").head(512)
with open("fna_paths.txt", "w") as f:
    for path in df["fna_path"]:
        f.write(path.strip() + "\n")
EOF

while IFS= read -r fna_gz; do
    if [[ -f "$fna_gz" ]]; then
        echo "Appending: $fna_gz"
        zcat "$fna_gz" >> "$COMBINED_FNA"
    else
        echo "Missing: $fna_gz"
    fi
done < fna_paths.txt

rm fna_paths.txt
echo "Combined .fna file created: $COMBINED_FNA"

# === STEP 2: Create BLAST DB ===
echo "Creating BLAST database..."
makeblastdb -in "$COMBINED_FNA" -dbtype nucl -out "$BLAST_DB"
echo "BLAST database created at: $BLAST_DB"
