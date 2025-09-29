import random
import os
import pandas as pd

# Define the base path
#dataset = "2"
#base_path = rf"/cs/usr/stavperez/sp/InSilicoSeq/dataset{dataset}"
base_path = rf"/cs/usr/stavperez/sp/InSilicoSeq/dataset_new"

# Load organism lists
#bact_df = pd.read_csv(f"ds{dataset}_shuffled_bact_org_info.csv")
#euk_df = pd.read_csv(f"ds{dataset}_shuffled_euk_org_info.csv")
bact_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/genome_refseq_bacteria.csv")
euk_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/genome_refseq_eukaryotes.csv")

# Select first mix_size organisms
mix_size = 128
selected_bact = bact_df["Organism"].head(mix_size).tolist()
selected_euk = euk_df["Organism"].head(mix_size).tolist()

def get_organism_paths(organism_list, org_type):
    """Generate file paths for selected organisms."""
    return {
        org: {
            "R1": f"{base_path}/output_reads_{org_type}/{org}_R1.fastq",
            "R2": f"{base_path}/output_reads_{org_type}/{org}_R2.fastq"
        } for org in organism_list
    }

# Create dictionaries for selected organisms
organism_data_bact = get_organism_paths(selected_bact, "Bact")
organism_data_euk = get_organism_paths(selected_euk, "Euk")

# Define output prefixes
output_prefix_all = f"{base_path}/output_reads_mix_256_orgs/output_reads_mix_256_orgs"
output_prefix_bact = f"{base_path}/output_reads_mix_256_orgs/output_reads_mix_128_bacts"
output_prefix_euk = f"{base_path}/output_reads_mix_256_orgs/output_reads_mix_128_euks"

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
    return paired_reads

def mix_multiple_organisms(organism_data, output_prefix):
    """Mix reads from multiple organisms into combined FASTQ files."""
    all_paired_reads = []
    print("Initial read counts:")
    for organism, files in organism_data.items():
        r1_count = count_reads(files['R1'])
        r2_count = count_reads(files['R2'])
#        print(f"{organism}: R1={r1_count}, R2={r2_count}")
        all_paired_reads.extend(process_organism_files(files['R1'], files['R2'], organism))

    random.shuffle(all_paired_reads)
    output_r1 = f"{output_prefix}_R1.fastq"
    output_r2 = f"{output_prefix}_R2.fastq"
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    with open(output_r1, 'w') as out_r1, open(output_r2, 'w') as out_r2:
        for r1, r2 in all_paired_reads:
            out_r1.write("\n".join(r1) + "\n")
            out_r2.write("\n".join(r2) + "\n")
    print(f"Output files: {output_r1}, {output_r2}")

# Mix reads and create final output files
mix_multiple_organisms({**organism_data_bact, **organism_data_euk}, output_prefix_all)
mix_multiple_organisms(organism_data_bact, output_prefix_bact)
mix_multiple_organisms(organism_data_euk, output_prefix_euk)
