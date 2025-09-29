#!/usr/bin/env python3
import pandas as pd
import os
import gzip
from Bio import SeqIO

# === Configuration ===
#base_dir = "/cs/usr/stavperez/sp/genomes_new/genomes_from_ncbi/"
base_dir = "/cs/usr/stavperez/sp/genomes_new/genomes_from_ncbi_unique_genus/"
#input_csv = "genomes_refseq_chromosome_unique_family.csv"
input_csv = "chromosome_assemblies_unique_genus_non_train.csv"
output_csv = "chromosome_assemblies_unique_genus_with_paths.csv"
bact_file = "genomes_bacteria_unique_genus.csv"
euk_file = "genomes_eukaryotes_unique_genus.csv"

# === Load Input CSV ===
df = pd.read_csv(input_csv)

# === Initialize Columns ===
fna_paths = []
gtf_paths = []
genome_lengths = []
sequence_headers = []
num_sequences = []

# === Process Each Genome ===
for gcf in df["assembly_accession"]:
    genome_dir = os.path.join(base_dir, gcf)
    fna_file = os.path.join(genome_dir, f"{gcf}_genomic.fna.gz")
    gtf_file = os.path.join(genome_dir, f"{gcf}_annotation.gtf.gz")

    fna_paths.append(fna_file if os.path.exists(fna_file) else None)
    gtf_paths.append(gtf_file if os.path.exists(gtf_file) else None)

    if os.path.exists(fna_file):
        try:
            total_len = 0
            headers = []
            with gzip.open(fna_file, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    total_len += len(record.seq)
                    headers.append(record.id)
            genome_lengths.append(total_len)
            sequence_headers.append(";".join(headers))
            num_sequences.append(len(headers))
        except Exception as e:
            print(f"[ERROR] Failed parsing {fna_file}: {e}")
            genome_lengths.append(None)
            sequence_headers.append(None)
            num_sequences.append(None)
    else:
        print(f"[WARN] Missing .fna.gz for {gcf}")
        genome_lengths.append(None)
        sequence_headers.append(None)
        num_sequences.append(None)

# === Add Computed Columns ===
df["fna_path"] = fna_paths
df["gtf_path"] = gtf_paths
df["genome_length"] = genome_lengths
df["sequence_headers"] = sequence_headers
df["num_sequences"] = num_sequences

# === Domain to Category Mapping ===
def map_domain(domain_val):
    if isinstance(domain_val, str):
        domain_val = domain_val.lower()
        if domain_val == "bacteria":
            return "bacteria"
        elif domain_val == "archaea":
            return "archaea"
        elif domain_val in {
            "fungi", "invertebrate", "plant", "protozoa",
            "vertebrate_mammalian", "vertebrate_other"
        }:
            return "eukaryotes"
    return "unknown"

df["category"] = df["domain"].apply(map_domain)

# === Save Full File ===
df.to_csv(output_csv, index=False)
print(f"✅ Full CSV saved to: {output_csv}")

# === Split by Category and Shuffle ===
#bacteria_df = df[df["category"] == "bacteria"].sample(frac=1, random_state=42).reset_index(drop=True)
#archaea_df = df[df["category"] == "archaea"].sample(frac=1, random_state=42).reset_index(drop=True)
#eukaryotes_df = df[df["category"] == "eukaryotes"].sample(frac=1, random_state=42).reset_index(drop=True)

# === Filter non-empty fna_path for category splits ===
bacteria_df = df[(df["category"] == "bacteria") & (df["fna_path"].notna())].sample(frac=1, random_state=42).reset_index(drop=True)
archaea_df = df[(df["category"] == "archaea") & (df["fna_path"].notna())].sample(frac=1, random_state=42).reset_index(drop=True)
eukaryotes_df = df[(df["category"] == "eukaryotes") & (df["fna_path"].notna())].sample(frac=1, random_state=42).reset_index(drop=True)


# === Save Shuffled Splits ===
bacteria_df.to_csv(bact_file, index=False)
archaea_df.to_csv("archaea.csv", index=False)
eukaryotes_df.to_csv(euk_file, index=False)

print("Shuffled CSVs created:")
print("  - genome_refseq_bacteria.csv")
print("  - genome_refseq_archaea.csv")
print("  - genome_refseq_eukaryotes.csv")
