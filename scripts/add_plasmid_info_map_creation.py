import pandas as pd
from Bio import SeqIO
import gzip
import os
from tqdm import tqdm

def extract_plasmid_labels(genomes_csv):
    df = pd.read_csv(genomes_csv)
    rows = []

    for _, row in tqdm(df.iterrows(), total=len(df)):
        organism = row["organism_name"]
        fna_path = row["fna_path"]

        if not os.path.exists(fna_path):
            print(f"[WARNING] Missing file: {fna_path}")
            continue

        open_func = gzip.open if fna_path.endswith(".gz") else open
        try:
            with open_func(fna_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    header = record.id  # or record.description for full header
                    is_plasmid = "plasmid" in record.description.lower()
                    rows.append({
                        "organism_name": organism,
                        "sequence_header": header,
                        "is_plasmid": is_plasmid
                    })
        except Exception as e:
            print(f"[ERROR] Failed reading {fna_path}: {e}")

    return pd.DataFrame(rows)

# Paths to your input files
bact_csv = "/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv"
euk_csv = "/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus.csv"

# Extract and combine
df_bact = extract_plasmid_labels(bact_csv)
df_euk = extract_plasmid_labels(euk_csv)
df_all = pd.concat([df_bact, df_euk], ignore_index=True)

# Save results
out_csv = "/cs/usr/stavperez/sp/genomes_new/sequence_header_to_plasmid.csv"
out_pkl = "/cs/usr/stavperez/sp/genomes_new/sequence_header_to_plasmid.pkl"
df_all.to_csv(out_csv, index=False)
df_all.to_pickle(out_pkl)

print(f"[DONE] Saved mapping for {len(df_all)} sequences")
