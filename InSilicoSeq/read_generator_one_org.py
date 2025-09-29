import os
import gzip
import shutil
import tempfile
import pandas as pd
from Bio import SeqIO
from iss import generator
import cProfile
import pstats
from io import StringIO
import argparse
import numpy as np
import time

def decompress_file(compressed_path, decompressed_path):
    """Decompress a gzip file to a new location."""
    print(f"Decompressing: {compressed_path} -> {decompressed_path}")
    with gzip.open(compressed_path, "rt") as fin, open(decompressed_path, "w") as fout:
        shutil.copyfileobj(fin, fout)

def get_genome_length(fasta_file):
    """
    Calculate the total genome length by reading the (possibly compressed) FASTA file.
    """
    print(f"Reading genome length from: {fasta_file}")
    total_length = 0
    # Open the file using gzip if it ends with .gz
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    with open_func(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_length += len(record.seq)
    print(f"Genome length: {total_length}")
    return total_length

def simulate_reads_for_organism(organism, fna_file, gff_file, output_dir, n_reads, err_mod, read_length=150):
    print(f"\n--- Simulating reads for {organism} ({n_reads:.0f} reads) ---")
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    # Prepare temporary decompressed files for fna and gff if necessary
    temp_fna = fna_file
    temp_gff = gff_file

    temp_files = []  # keep track of temporary files to delete later

    # Use the process ID to ensure unique temporary file names
    pid = os.getpid()

    if fna_file.endswith('.gz'):
        temp_fna = os.path.join(tempfile.gettempdir(), f"{organism}_{pid}_temp.fasta")
        decompress_file(fna_file, temp_fna)
        temp_files.append(temp_fna)

    if gff_file and gff_file.endswith('.gz'):
        temp_gff = os.path.join(tempfile.gettempdir(), f"{organism}_{pid}_temp.gff")
        decompress_file(gff_file, temp_gff)
        temp_files.append(temp_gff)

    forward_handle = open(os.path.join(output_dir, f"{organism}_R1.fastq"), "w")
    reverse_handle = open(os.path.join(output_dir, f"{organism}_R2.fastq"), "w")
    mutations_handle = open(os.path.join(output_dir, f"{organism}.tsv"), "w")

    try:
        print(f"Calling 'simulate_and_add_info_to_reads' for {organism}")

        generator.simulate_and_add_info_to_reads(
            temp_gff, temp_fna, err_mod,
            n_pairs=int(n_reads), cpu_number=4,
            forward_handle=forward_handle,
            reverse_handle=reverse_handle,
            mutations_handle=mutations_handle,
            sequence_type="metagenomics",
            gc_bias=True
        )
        print(f"Read simulation completed successfully for {organism}")
    except Exception as e:
        print(f"ERROR: Simulation failed for {organism}: {e}")
    finally:
        forward_handle.close()
        reverse_handle.close()
        mutations_handle.close()
        # Clean up temporary files if they were created
        for tmp in temp_files:
            if os.path.exists(tmp):
                os.remove(tmp)

    end_time = time.time()
    print(f"Read simulation for {organism} completed in {end_time - start_time:.2f} seconds")


def generate_reads_from_csv(csv_file, org_type, n_reads, err_mod,
                            read_length=150, output_base="."):

    print(f"\n--- Loading CSV file: {csv_file} ---")
    if not os.path.exists(csv_file):
        print(f"ERROR: CSV file '{csv_file}' not found.")
        exit(1)
    df = pd.read_csv(csv_file)
    print(f"Loaded CSV with {len(df)} rows.")

    #output_dir = f"output_reads_{org_type}"
    output_dir = os.path.join(output_base, f"output_reads_{org_type}")
    os.makedirs(output_dir, exist_ok=True)


    for _, row in df.iterrows():
        organism = row["organism_name"]
        fna_file = row["fna_path"]
        gff_file = row["gtf_path"]
        
        print(f"\nProcessing {organism} - FNA: {fna_file}, GFF: {gff_file}")

        if not os.path.exists(fna_file):
            print(f"ERROR: FNA file {fna_file} does not exist!")
            continue

        if gff_file and not os.path.exists(gff_file):
            print(f"WARNING: GFF file {gff_file} does not exist! Skipping GFF.")


      #  if org_type == "Bact":
     #       genome_length = row["genome_length"]
    #        n_reads = (genome_length)*3 / 300
   #         simulate_reads_for_organism(organism, fna_file, gff_file, output_dir, n_reads, err_mod, read_length)

    # Process bacteria genomes if available
  #      if org_type == "Euk":
#            simulate_reads_for_organism(organism, fna_file, gff_file,
 #                                       output_dir, n_reads, err_mod, read_length)
        simulate_reads_for_organism(organism, fna_file, gff_file,
                                        output_dir, n_reads, err_mod,
                                        read_length)



def profile_generate_reads(euk_genomes, bact_genomes, err_mod):
    pr = cProfile.Profile()
    pr.enable()
    generate_reads_from_csv(euk_genomes, bact_genomes, err_mod)
    pr.disable()

    s = StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
    ps.print_stats(20)  # Print the top 20 functions
    print(s.getvalue())


def create_organisms_dict(org="Euk"):
    # Base directory containing all the RefSeq genomes
    base_dir = f"/sci/labs/aerez/stavperez/genomes_{org}_from_ncbi/refseq"
    org_genomes = {}
    org_annotations = {}
    counter = 1

    # Traverse each category directory (e.g., fungi, invertebrate, plant, etc.)
    for category in os.listdir(base_dir):
        category_path = os.path.join(base_dir, category)
        if os.path.isdir(category_path):
            # Loop over each genome directory inside the category
            for genome_dir in os.listdir(category_path):
                genome_path = os.path.join(category_path, genome_dir)
                if os.path.isdir(genome_path):
                    fna_file = None
                    gff_file = None
                    # Look for the compressed fna and gff files in the genome directory
                    for file in os.listdir(genome_path):
                        if file.endswith("_genomic.fna.gz"):
                            fna_file = os.path.join(genome_path, file)
                        elif file.endswith("_genomic.gff.gz"):
                            gff_file = os.path.join(genome_path, file)
                    # Only add to the dictionaries if the fna file exists
                    if fna_file:
                        key = f"{org}_{counter}"
                        org_genomes[key] = {"fna": fna_file, "gff": gff_file}
                        org_annotations[key] = genome_dir
                        counter += 1
                        # Break immediately if we've reached 640 organisms
                        if counter > 640:
                            return org_genomes, org_annotations
    return org_genomes, org_annotations


def write_org_annotations_excel(org_annotations, output_file):
    """
    Create an Excel file with two columns: Organism and annotation,
    containing the contents of the org_annotations dictionary.
    """
    df = pd.DataFrame(list(org_annotations.items()), columns=["Organism", "annotation"])
    df.to_csv(output_file, index=False)
    print(f"Excel file created: {output_file}")

def write_org_info_csv(org_genomes, org_annotations, output_file):
    """
    Create a CSV file that contains:
      - The organism key (e.g. "Euk_1")
      - Its annotation (from org_annotations)
      - The fna and gff file paths (from euk_genomes)
      - The genome length computed from the fna file using get_genome_length.
    """
    data = []
    for org_id, genome_info in org_genomes.items():
        # Get the corresponding annotation (if available)
        annotation = org_annotations.get(org_id, "")
        fna = genome_info.get("fna")
        gff = genome_info.get("gff")
        # Compute genome length using get_genome_length
        genome_length = get_genome_length(fna) if fna else None
        data.append({
            "Organism": org_id,
            "annotation": annotation,
            "fna": fna,
            "gff": gff,
            "genome_length": genome_length
        })
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"CSV file created: {output_file}")





if __name__ == "__main__":
    # org = "Bact"
    # org_genomes, org_annotations = create_organisms_dict(org)
    # write_org_annotations_excel(org_annotations, f"{org}_annotations.csv")
    # write_org_info_csv(org_genomes, org_annotations, f"{org}_org_info.csv")
    #
    parser = argparse.ArgumentParser(description="Run read simulation for a specific organism from CSV.")
    parser.add_argument("--csv", type=str, required=True, help="Path to the CSV file containing organism information.")
    parser.add_argument("--org", type=str, required=True, choices=["Euk", "Bact"],
                    help="Organism type: Euk or Bact")
    parser.add_argument("--output_dir", type=str, required=True,
                    help="Path to the base output directory")


    args = parser.parse_args()
    csv_file = args.csv
    print(f"\n### Starting Read Generation for CSV: {csv_file} ###")

    if not os.path.exists(csv_file):
        print(f"Error: CSV file '{csv_file}' not found.")
        exit(1)

    print(f"Loading CSV file: {csv_file}")

    print("Loading error model...")
    err_mod = generator.load_error_model("kde", 0, "miseq", None, None, True)
    print("Error model loaded successfully.")
    #    generate_reads(euk_genomes, bact_genomes, err_mod)
#     profile_generate_reads("euk_org_info.csv", bact_genomes, err_mod)
    # generate_reads_from_csv("euk_org_info.csv", bact_genomes, err_mod)
    # n_reads = 12068158
    n_reads = 35000
#    n_reads = 100000

    generate_reads_from_csv(csv_file, args.org, n_reads, err_mod, read_length=150, output_base=args.output_dir)

    #generate_reads_from_csv(csv_file, "Euk", n_reads, err_mod, read_length=150)
#    generate_reads_from_csv(csv_file, "Bact", n_reads, err_mod, read_length=150)

    print("\n### Read Generation Process Completed ###")


    # simulate_reads_for_organism("Bact_1", fna_file, gff_file, output_dir,
    #                             n_reads, err_mod, read_length=150)
