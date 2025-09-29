import random
import os
import pandas as pd

# Base path where the reads are stored
base_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus"

# Load organism lists
bact_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv")
euk_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus.csv")

# Select top 128*3 from each and split into 3 chunks of 128
print("Loading organism lists...")
bact_chunks = [bact_df["organism_name"].iloc[i*128:(i+1)*128].tolist() for i in range(3)]
euk_chunks = [euk_df["organism_name"].iloc[i*128:(i+1)*128].tolist() for i in range(3)]

###################
print("\n=== Verifying Euk chunk 3 content ===")
expected_rows = euk_df["organism_name"].iloc[256:384].tolist()
actual_chunk = euk_chunks[2]

# Print a few from each to compare
print("Expected (from iloc 256:384):")
print(expected_rows[:5], "...", expected_rows[-5:])

print("Actual chunk 3 (euk_chunks[2]):")
print(actual_chunk[:5], "...", actual_chunk[-5:])

# Check if all match
if expected_rows == actual_chunk:
    print("✅ Chunk 3 correctly contains rows 256–383 from the CSV.")
else:
    print("❌ Mismatch detected in chunk 3! Something is wrong.")
################

print(f"Bacterial organisms loaded: {len(bact_df)}, Eukaryotic organisms loaded: {len(euk_df)}")


def get_organism_paths(organism_list, org_type):
    print(organism_list)
    """Generate file paths for selected organisms."""
    print(f"Generating paths for {len(organism_list)} {org_type} organisms...")
    return {
        org: {
            "R1": f"{base_path}/output_reads_{org_type}/{org}_R1.fastq",
            "R2": f"{base_path}/output_reads_{org_type}/{org}_R2.fastq"
        } for org in organism_list
    }

def count_reads(file_path: str) -> int:
    """Count the number of reads in a FASTQ file."""
    try:
        with open(file_path, 'r') as f:
            return sum(1 for _ in f) // 4
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return 0

def process_organism_files(r1_file: str, r2_file: str, organism_name: str):
    """Process paired-end FASTQ files for a single organism."""
    paired_reads = []
    try:
        with open(r1_file, 'r') as r1, open(r2_file, 'r') as r2:
            while True:
                r1_lines = [r1.readline().strip() for _ in range(4)]
                r2_lines = [r2.readline().strip() for _ in range(4)]

                if not r1_lines[0] or not r2_lines[0]:
                    break

                r1_lines[0] += f" Source:{organism_name}"
                r2_lines[0] += f" Source:{organism_name}"

                paired_reads.append((r1_lines, r2_lines))
    except FileNotFoundError:
        print(f"Missing files for {organism_name}: {r1_file} or {r2_file}")
    print(f"{organism_name}: {len(paired_reads)} paired reads processed.")
    return paired_reads

def mix_multiple_organisms(organism_data, output_prefix):
    """Mix reads from multiple organisms into combined FASTQ files."""
    all_paired_reads = []
    print(f"Processing {len(organism_data)} organisms...")
    for organism, files in organism_data.items():
        all_paired_reads.extend(process_organism_files(files['R1'], files['R2'], organism))

    random.shuffle(all_paired_reads)
    output_r1 = f"{output_prefix}_R1.fastq"
    output_r2 = f"{output_prefix}_R2.fastq"
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    with open(output_r1, 'w') as out_r1, open(output_r2, 'w') as out_r2:
        for r1, r2 in all_paired_reads:
            out_r1.write("\n".join(r1) + "\n")
            out_r2.write("\n".join(r2) + "\n")
    print(f"Output completed: {output_r1}, {output_r2}, Total pairs: {len(all_paired_reads)}")


# Process and mix 3 groups
#for i in range(3):
 #   if i == 0 or i == 1:
  #      continue
   # bact_group = bact_chunks[i]
    #euk_group = euk_chunks[i]
 #   output_prefix = f"{base_path}/output_reads_mix_256_orgs/mix_{i+1}_128bact_128euk"
  #  organism_data = {**get_organism_paths(bact_group, "Bact"), **get_organism_paths(euk_group, "Euk")}
   # mix_multiple_organisms(organism_data, output_prefix)


for i in range(3):
    if i == 0 or i == 1:
        continue

    print(f"\n=== Preparing dataset {i+1} ===")
    bact_group = bact_chunks[i]
    euk_group = euk_chunks[i]

    print(f"  Bacteria chunk {i+1} ({len(bact_group)}):")
    for org in bact_group:
        print(f"    Bact: {org} from chunk {i+1}")

    print(f"  Eukaryote chunk {i+1} ({len(euk_group)}):")
    for org in euk_group:
        print(f"    Euk: {org} from chunk {i+1}")

    output_prefix = f"{base_path}/output_reads_mix_256_orgs/mix_{i+1}_128bact_128euk"

    organism_data = {
        **get_organism_paths(bact_group, "Bact"),
        **get_organism_paths(euk_group, "Euk")
    }

#    mix_multiple_organisms(organism_data, output_prefix)
