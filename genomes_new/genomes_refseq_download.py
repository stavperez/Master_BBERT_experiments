#!/usr/bin/env python3
import os
import csv
import requests
from Bio import Entrez

Entrez.email = "stav.perez@mail.huji.ac.il"

def get_ftp_path(refseq_id):
    """Get FTP path from RefSeq GCF ID using Entrez."""
    try:
        search = Entrez.esearch(db="assembly", term=refseq_id, retmode="xml")
        record = Entrez.read(search)
        if not record["IdList"]:
            print(f"[WARN] No Entrez record found for {refseq_id}")
            return None
        uid = record["IdList"][0]

        summary_handle = Entrez.esummary(db="assembly", id=uid, retmode="xml")
        summary_data = Entrez.read(summary_handle)
        docsum = summary_data["DocumentSummarySet"]["DocumentSummary"][0]

        ftp_path = docsum.get("FtpPath_RefSeq")
        return ftp_path

    except Exception as e:
        print(f"[ERROR] Failed to get FTP path for {refseq_id}: {e}")
        return None

def download_file(url, out_path):
    """Download a file from a URL and save it."""
    url = url.replace("ftp://", "https://")
    try:
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(1024):
                    f.write(chunk)
            print(f"[OK] Downloaded {os.path.basename(out_path)}")
        else:
            print(f"[WARN] Failed to download {url} (status {r.status_code})")
    except Exception as e:
        print(f"[ERROR] Exception downloading {url}: {e}")

def download_by_refseq_id(refseq_id, output_dir):
    ftp_path = get_ftp_path(refseq_id)
    if not ftp_path:
        print(f"[SKIP] No FTP path found for {refseq_id}")
        return

    asm_name = os.path.basename(ftp_path)
    base_url = ftp_path.replace("ftp://", "https://")
    files = {
        "fna": f"{base_url}/{asm_name}_genomic.fna.gz",
        "gff": f"{base_url}/{asm_name}_genomic.gff.gz",
        "gtf": f"{base_url}/{asm_name}_genomic.gtf.gz"
    }

    org_out_dir = os.path.join(output_dir, refseq_id)
    os.makedirs(org_out_dir, exist_ok=True)

    # Download genome
    download_file(files["fna"], os.path.join(org_out_dir, f"{refseq_id}_genomic.fna.gz"))

    # Download annotation (prefer GTF, fallback to GFF)
    gtf_path = os.path.join(org_out_dir, f"{refseq_id}_annotation.gtf.gz")
    gff_path = os.path.join(org_out_dir, f"{refseq_id}_annotation.gff.gz")
    r = requests.get(files["gtf"], stream=True)
    if r.status_code == 200:
        with open(gtf_path, "wb") as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)
        print(f"[OK] Downloaded GTF annotation for {refseq_id}")
    else:
        download_file(files["gff"], gff_path)

def main():
#    csv_path = "genomes_refseq_chromosome_unique_family.csv"
    csv_path = "chromosome_assemblies_unique_genus_non_train.csv"
    output_dir = "genomes_from_ncbi_unique_genus"
    os.makedirs(output_dir, exist_ok=True)

    bacteria_count = 0
    max_bacteria = 752

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            domain = row.get("domain", "").strip().lower()
            refseq_id = row.get("assembly_accession", "").strip()
            
            # Skip archaea
            if domain == "archaea":
                continue

            # If bacteria, enforce 752 limit
            if domain == "bacteria":
                if bacteria_count >= max_bacteria:
                    continue
                bacteria_count += 1

            # Process other organisms or allowed bacteria
            if refseq_id:
                print(f"\n--- Processing {refseq_id} ---")
                download_by_refseq_id(refseq_id, output_dir)

if __name__ == "__main__":
    main()
