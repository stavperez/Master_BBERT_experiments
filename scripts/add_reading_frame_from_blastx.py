import os
import subprocess
import pandas as pd
from Bio import SeqIO

# === Paths ===
fastq_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_Bact/Pseudomonas_aeruginosa_R1.fastq"
#csv_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_aeruginosa/Pseudomonas_aeruginosa_go_terms.csv"
csv_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_Saccharomyces_GO_info.csv"
faa_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_aeruginosa/protein.faa"

# === Temporary working directory ===
tmp_dir = "tmp_blastx"
os.makedirs(tmp_dir, exist_ok=True)

fasta_path = os.path.join(tmp_dir, "reads.fasta")
blast_db_prefix = os.path.join(tmp_dir, "blastx_db")
blast_output_path = os.path.join(tmp_dir, "blastx_output.tsv")
output_csv_path = csv_path.replace(".csv", "_with_frame.csv")

# === Step 1: Convert FASTQ to FASTA ===
print(f"🔄 Converting FASTQ to FASTA: {fastq_path}")
try:
    count = 0
    with open(fasta_path, "w") as out_fasta:
        for record in SeqIO.parse(fastq_path, "fastq"):
            SeqIO.write(record, out_fasta, "fasta")
            count += 1
    print(f"✅ FASTA file created: {fasta_path} with {count} reads")
except Exception as e:
    print(f"❌ Error during FASTQ to FASTA conversion: {e}")
    exit(1)

# === Step 2: Create BLAST protein DB ===
print(f"🔬 Creating BLAST database from: {faa_path}")
if not os.path.exists(faa_path):
    print("❌ protein.faa file not found. Aborting.")
    exit(1)

try:
    subprocess.run([
        "makeblastdb",
        "-in", faa_path,
        "-dbtype", "prot",
        "-out", blast_db_prefix
    ], check=True)
    print(f"✅ BLAST DB created at: {blast_db_prefix}")
except subprocess.CalledProcessError as e:
    print(f"❌ Failed to create BLAST DB: {e}")
    exit(1)

# === Step 3: Run BLASTX with sframe ===
print("🚀 Running BLASTX...")
try:
    blastx_cmd = [
    "blastx",
    "-query", fasta_path,
    "-db", blast_db_prefix,
    "-evalue", "1e-5",
    "-max_target_seqs", "1",
    "-out", blast_output_path,
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe"
]
    subprocess.run(blastx_cmd, check=True)
    print(f"✅ BLASTX finished. Output saved to: {blast_output_path}")
except subprocess.CalledProcessError as e:
    print(f"❌ BLASTX failed: {e}")
    exit(1)

# === Step 4: Load BLAST results ===
print(f"📥 Loading BLAST results from: {blast_output_path}")
if not os.path.exists(blast_output_path) or os.path.getsize(blast_output_path) == 0:
    print("❌ No BLASTX output found. Something went wrong.")
    exit(1)

blast_cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "reading_frame"
]

try:
    blast_df = pd.read_csv(blast_output_path, sep="\t", names=blast_cols)
    print(f"✅ Loaded BLASTX hits: {len(blast_df)} rows")
    print("🔍 BLAST sample:")
    print(blast_df.head(3))
    print(f"🧪 Unique read IDs in BLAST: {blast_df['qseqid'].nunique()}")
except Exception as e:
    print(f"❌ Failed to load BLAST output: {e}")
    exit(1)

# === Step 5: Filter best hit per read ===
print("🧠 Selecting best hit per read...")
blast_best = (
    blast_df.sort_values("bitscore", ascending=False)
    .drop_duplicates("qseqid")[["qseqid", "reading_frame"]]
    .rename(columns={"qseqid": "read_id", "reading_frame": "reading_frame"})
)
print(f"✅ Best reading frame selected for {len(blast_best)} reads")
print("🔍 Best BLAST sample:")
print(blast_best.head(3))

# === Step 6: Load and update original CSV ===
print(f"📄 Loading original CSV: {csv_path}")
try:
    df = pd.read_csv(csv_path)
    print(f"✅ Loaded: {len(df)} rows")
    print("🔍 Sample read_ids from CSV:")
    print(df["read_id"].drop_duplicates().head(3).to_list())
except Exception as e:
    print(f"❌ Failed to load CSV: {e}")
    exit(1)

# === Normalize read IDs ===
print("🛠️ Normalizing read IDs for merge...")
blast_best["read_id_normalized"] = blast_best["read_id"].astype(str).str.replace(r'/[12]$', '', regex=True)
df["read_id_normalized"] = df["read_id"].astype(str).str.replace(r'/[12]$', '', regex=True)

# === Merge
merged = df.merge(blast_best[["read_id_normalized", "reading_frame"]],
                  on="read_id_normalized", how="left")

merged.drop(columns=["read_id_normalized"], inplace=True)

# === Final stats
found = merged["reading_frame"].notna().sum()
missing = merged["reading_frame"].isna().sum()
print(f"✅ Frame added to {found} reads, {missing} reads missing (no hit)")

print("\n🔍 Sample merged rows with frame:")
print(merged[merged["reading_frame"].notna()].head(5))

# === Step 7: Save output CSV ===
print(f"💾 Saving updated CSV to: {output_csv_path}")
try:
    merged.to_csv(output_csv_path, index=False)
    print("✅ CSV saved successfully.")
except Exception as e:
    print(f"❌ Failed to save CSV: {e}")
