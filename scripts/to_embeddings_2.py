"""
This script generates balanced FASTQ datasets from a labeled CSV file.

1. For reading frames +1, -1, and +2:
   - Selects an equal number of TP and FN reads (50/50 split).
   - Saves matching R1 and R2 FASTQ files per frame.

2. For coding/non-coding reads:
   - Selects equal number of 'coding' and 'non-coding' reads based on the "CDS" column.
   - Saves matching R1 and R2 FASTQ files.

Output FASTQs are saved to a specified directory.
"""

import os
import pandas as pd
from Bio import SeqIO

# === INPUT PATHS ===
fastq_path_R1 = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact/Pseudomonas_aeruginosa_R1.fastq"
fastq_path_R2 = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact/Pseudomonas_aeruginosa_R2.fastq"
csv_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_Saccharomyces_GO_info_with_frame.csv"

# === OUTPUT DIRECTORY ===
output_dir = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/balanced_fastqs"
os.makedirs(output_dir, exist_ok=True)

# === PARAMETERS ===
frames_of_interest = [1, -1, 2]
total_reads_per_set = 1000  # 500 TP + 500 FN, 500 coding + 500 non-coding

# === Load CSV once ===
import pandas as pd

# === Input CSV path ===
csv_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_Saccharomyces_GO_info_with_frame.csv"
df = pd.read_csv(csv_path)

# === 1. Total counts ===
print("=== TOTAL COUNTS ===")
print(f"Total reads:         {len(df)}")
print(f"Total TP:            {(df['confusion_type'] == 'TP').sum()}")
print(f"Total FN:            {(df['confusion_type'] == 'FN').sum()}")
print(f"Total CODING:        {(df['CDS'] == 'coding').sum()}")
print(f"Total NON-CODING:    {(df['CDS'] == 'non-coding').sum()}")

# === 2. TP and FN per reading frame (including None) ===
print("\n=== TP per reading frame (including None) ===")
tp_rf = df[df["confusion_type"] == "TP"]["reading_frame"].fillna("None").astype(str)
print(tp_rf.value_counts().sort_index())

print("\n=== FN per reading frame (including None) ===")
fn_rf = df[df["confusion_type"] == "FN"]["reading_frame"].fillna("None").astype(str)
print(fn_rf.value_counts().sort_index())

# === 3. CODING / NON-CODING per reading frame (including None) ===
print("\n=== CODING / NON-CODING per reading frame (including None) ===")
df["reading_frame_str"] = df["reading_frame"].fillna("None").astype(str)
coding_stats = df.groupby(["reading_frame_str", "CDS"]).size().unstack(fill_value=0)
print(coding_stats)



# === Task 1: For each frame, generate 50/50 TP/FN FASTQ ===
for frame in frames_of_interest:
    df_frame = df[df["reading_frame"] == frame]
    tp_ids = df_frame[df_frame["confusion_type"] == "TP"]["read_id"].tolist()
    fn_ids = df_frame[df_frame["confusion_type"] == "FN"]["read_id"].tolist()

    min_len = min(len(tp_ids), len(fn_ids), total_reads_per_set // 2)
    tp_ids_set = set(tp_ids[:min_len])
    fn_ids_set = set(fn_ids[:min_len])
    selected_ids = tp_ids_set.union(fn_ids_set)
    selected_ids_R2 = {rid.replace("/1", "/2") for rid in selected_ids}

    print(f"✅ Frame {frame}: {min_len} TP + {min_len} FN reads selected")

    # === Extract R1 ===
    frame_R1_reads = []
    with open(fastq_path_R1, "r") as f1:
        for record in SeqIO.parse(f1, "fastq"):
            if record.id in selected_ids:
                frame_R1_reads.append(record)

    # === Extract R2 ===
    frame_R2_reads = []
    with open(fastq_path_R2, "r") as f2:
        for record in SeqIO.parse(f2, "fastq"):
            if record.id in selected_ids_R2:
                frame_R2_reads.append(record)

    # === Save FASTQs ===
    out_r1 = os.path.join(output_dir, f"reading_frame_{frame}_balanced_confusion_R1.fastq")
    out_r2 = os.path.join(output_dir, f"reading_frame_{frame}_balanced_confusion_R2.fastq")
    SeqIO.write(frame_R1_reads, out_r1, "fastq")
    SeqIO.write(frame_R2_reads, out_r2, "fastq")
    print(f"   → Saved to {out_r1} and {out_r2}")

# === Task 2: Balanced coding/non-coding ===
coding_ids = df[df["CDS"] == "coding"]["read_id"].tolist()
noncoding_ids = df[df["CDS"] == "non-coding"]["read_id"].tolist()
min_len_coding = min(len(coding_ids), len(noncoding_ids), total_reads_per_set // 2)
coding_ids_set = set(coding_ids[:min_len_coding])
noncoding_ids_set = set(noncoding_ids[:min_len_coding])
cds_ids = coding_ids_set.union(noncoding_ids_set)
cds_ids_R2 = {rid.replace("/1", "/2") for rid in cds_ids}

# === Extract R1 ===
coding_R1_reads = []
with open(fastq_path_R1, "r") as f1:
    for record in SeqIO.parse(f1, "fastq"):
        if record.id in cds_ids:
            coding_R1_reads.append(record)

# === Extract R2 ===
coding_R2_reads = []
with open(fastq_path_R2, "r") as f2:
    for record in SeqIO.parse(f2, "fastq"):
        if record.id in cds_ids_R2:
            coding_R2_reads.append(record)

# === Save ===
out_coding_r1 = os.path.join(output_dir, "balanced_coding_R1.fastq")
out_coding_r2 = os.path.join(output_dir, "balanced_coding_R2.fastq")
SeqIO.write(coding_R1_reads, out_coding_r1, "fastq")
SeqIO.write(coding_R2_reads, out_coding_r2, "fastq")
print(f"✅ Coding/non-coding: Saved {len(coding_R1_reads)} R1 reads and {len(coding_R2_reads)} R2 reads")
