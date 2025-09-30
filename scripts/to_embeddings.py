"""
This script filters reads from a FASTQ file and CSV labels to generate:

1. FASTQ files grouped by reading frame (-3 to +3).
2. A combined FASTQ with an equal number of TP and FN reads (50/50 split).
3. Matching R2 FASTQ files for all selected R1 reads.

Used for organizing input data for downstream embedding or classification.
"""


import os
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# === INPUT FILES ===
fastq_path_R1 = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact/Pseudomonas_aeruginosa_R1.fastq"
fastq_path_R2 = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact/Pseudomonas_aeruginosa_R2.fastq"
csv_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_Saccharomyces_GO_info_with_frame.csv"

# === OUTPUT FILES ===
output_dir_frames = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/reads_by_reading_frame/"
output_confusion_R1_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/reads_by_confusion_type.fastq"
output_confusion_R2_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/reads_by_confusion_type_R2.fastq"

# === PARAMETERS ===
frame_limit = 2048
confusion_limit = 1024
valid_frames = {-3, -2, -1, 1, 2, 3}

# === Setup output directory ===
os.makedirs(output_dir_frames, exist_ok=True)

# === Load CSV ===
df = pd.read_csv(csv_path)

# === Create mappings ===
frame_map = df[df["reading_frame"].isin(valid_frames)].set_index("read_id")["reading_frame"].to_dict()
confusion_TP = set(df[df["confusion_type"] == "TP"]["read_id"])
confusion_FN = set(df[df["confusion_type"] == "FN"]["read_id"])

# === Prepare containers ===
frame_selected = defaultdict(list)
frame_selected_ids = defaultdict(set)
confusion_selected = []
confusion_selected_ids = set()

# === Parse R1 FASTQ and collect reads ===
with open(fastq_path_R1, "r") as f1:
    for record in SeqIO.parse(f1, "fastq"):
        read_id = record.id

        # For reading frame groups
        if read_id in frame_map:
            frame = frame_map[read_id]
            if len(frame_selected[frame]) < frame_limit:
                frame_selected[frame].append(record)
                frame_selected_ids[frame].add(read_id)

        # For confusion type groups
        if read_id in confusion_TP and sum(1 for r in confusion_selected if r[1] == "TP") < confusion_limit:
            confusion_selected.append((record, "TP"))
            confusion_selected_ids.add(read_id)

        elif read_id in confusion_FN and sum(1 for r in confusion_selected if r[1] == "FN") < confusion_limit:
            confusion_selected.append((record, "FN"))
            confusion_selected_ids.add(read_id)

        if all(len(frame_selected[fr]) >= frame_limit for fr in valid_frames) and \
           sum(1 for r in confusion_selected if r[1] == "TP") >= confusion_limit and \
           sum(1 for r in confusion_selected if r[1] == "FN") >= confusion_limit:
            break

# === Write R1 FASTQ files ===
for frame, records in frame_selected.items():
    out_path = os.path.join(output_dir_frames, f"reading_frame_{frame}_R1.fastq")
    with open(out_path, "w") as f_out:
        SeqIO.write(records, f_out, "fastq")
    print(f"✅ Saved {len(records)} reads to {out_path}")

with open(output_confusion_R1_path, "w") as f_out:
    SeqIO.write([rec for rec, _ in confusion_selected], f_out, "fastq")
print(f"✅ Saved {len(confusion_selected)} reads to {output_confusion_R1_path}")

# === Prepare target R2 read IDs ===
frame_selected_ids_R2 = {frame: {rid.replace("/1", "/2") for rid in rids}
                         for frame, rids in frame_selected_ids.items()}
confusion_selected_ids_R2 = {rid.replace("/1", "/2") for rid in confusion_selected_ids}

# === Parse R2 FASTQ and write matching reads ===
frame_output_R2_handles = {
    frame: open(os.path.join(output_dir_frames, f"reading_frame_{frame}_R2.fastq"), "w")
    for frame in valid_frames
}
confusion_output_R2 = open(output_confusion_R2_path, "w")

with open(fastq_path_R2, "r") as f2:
    for record in SeqIO.parse(f2, "fastq"):
        rid = record.id

        # Check reading frame groups
        for frame in valid_frames:
            if rid in frame_selected_ids_R2[frame]:
                SeqIO.write(record, frame_output_R2_handles[frame], "fastq")

        # Check confusion
        if rid in confusion_selected_ids_R2:
            SeqIO.write(record, confusion_output_R2, "fastq")

# === Close R2 output files ===
for f in frame_output_R2_handles.values():
    f.close()
confusion_output_R2.close()

print(f"✅ Saved all corresponding R2 reads.")
