import os
import pandas as pd
import requests
from Bio import Entrez
from io import StringIO
from time import sleep
from tqdm import tqdm

# === Configuration ===
Entrez.email = "stav.perez@mail.huji.ac.il"
REFSEQ_FTPS = {
    "bacteria": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "archaea": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
    "fungi": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
    "protozoa": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt",
    "plant": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt",
    "invertebrate": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/assembly_summary.txt",
    "vertebrate_mammalian": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt",
    "vertebrate_other": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/assembly_summary.txt"
}

TAXONOMIC_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']

# === Functions ===
def download_and_parse_summary(name, url):
    print(f"Downloading {name}...")
    response = requests.get(url)
    response.raise_for_status()

    lines = response.text.splitlines()
    header_line = [line for line in lines if line.startswith('#')][-1]
    header = header_line.lstrip('# ').split('\t')
    data_lines = [line for line in lines if not line.startswith('#')]

    df = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t', names=header)
    df["domain"] = name

    required_cols = ['assembly_accession', 'organism_name', 'taxid', 'assembly_level']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing columns in {name}: {df.columns}")

    return df[required_cols + ['domain']]


def get_lineage(taxid):
    try:
        handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
        records = Entrez.read(handle)
        lineage_info = {rank: None for rank in TAXONOMIC_RANKS}
        if records:
            lineage = records[0].get("LineageEx", [])
            for rank in lineage:
                if rank['Rank'] in lineage_info:
                    lineage_info[rank['Rank']] = rank['ScientificName']
            lineage_info['species'] = records[0].get("ScientificName", None)
        return lineage_info
    except Exception as e:
        print(f"Error for taxid {taxid}: {e}")
        return {rank: None for rank in TAXONOMIC_RANKS}

# === Main ===
def main():
    all_dfs = []
    for name, url in REFSEQ_FTPS.items():
        try:
            df = download_and_parse_summary(name, url)
            all_dfs.append(df)
        except Exception as e:
            print(f"Failed to process {name}: {e}")

    print("Merging all domains...")
    df_all = pd.concat(all_dfs, ignore_index=True).drop_duplicates()

    # Save full list WITHOUT taxonomy
    df_all.to_csv("all_assemblies.csv", index=False)
    print("✅ Saved: all_assemblies.csv (without taxonomy)")

    # Filter Chromosome + Complete Genome
    keep_levels = ["chromosome", "complete genome"]
    df_chr = df_all[df_all['assembly_level'].str.lower().isin(keep_levels)].copy()

    print(f"Fetching taxonomic lineages for {len(df_chr)} chromosome/complete genome assemblies...")
    unique_taxids = df_chr['taxid'].dropna().astype(int).unique()

    taxonomy_records = []
    for taxid in tqdm(unique_taxids):
        lineage = get_lineage(taxid)
        lineage['taxid'] = taxid
        taxonomy_records.append(lineage)
        sleep(0.3)  # be kind to NCBI

    df_tax = pd.DataFrame(taxonomy_records)
    df_chr_merged = df_chr.merge(df_tax, on='taxid', how='left')

    df_chr_merged.to_csv("chromosome_assemblies_with_taxonomy.csv", index=False)
    print("✅ Saved: chromosome_assemblies_with_taxonomy.csv (with taxonomy)")


if __name__ == "__main__":
    main()
