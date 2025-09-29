#!/usr/bin/env python3
import pandas as pd
import subprocess
import re
import os

def check_genome_availability(organism_name):
    """
    Uses ncbi-genome-download in dry-run mode to check if a RefSeq genome 
    for the given organism is available in the desired formats (.fna and .gtf).
    """
    result = subprocess.run([
        "ncbi-genome-download", "--dry-run", "--formats", "fna,gtf",
        "all", "-s", "refseq", "-g", organism_name
    ], capture_output=True, text=True)
    # Search the dry-run output for a valid genome ID (e.g. GCF_000001405.39)
    print(f"Results for {organism_name}:")
    print(result.stdout)
    if re.search(r'GCF_\d+\.\d+', result.stdout):
        return True
    else:
        return False

def download_genome(organism_name, output_dir="genomes_proc_from_ncbi"):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    command = [
        "ncbi-genome-download", "all", "-s", "refseq",
        "--formats", "fasta,gff", "-g", organism_name,
        "-o", output_dir
    ]
    subprocess.run(command)
    print(f"Downloaded genome for {organism_name}")

def main():
    # Update this with your CSV file path; ensure that the CSV has a column called 'Organism'
    csv_file = "/cs/usr/stavperez/sp/prokaryotes_512_organisms_with_genomes.csv"
    df = pd.read_csv(csv_file)
    
    for organism in df["Organism Name"]:
        print(f"Processing {organism}...")
#        if check_genome_availability(organism):
        download_genome(organism)
#        else:
 #           print(f"Genome for {organism} not available in RefSeq.")

if __name__ == "__main__":
    main()
