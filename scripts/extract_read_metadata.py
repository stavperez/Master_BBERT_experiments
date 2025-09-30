# extract_read_metadata.py
import pandas as pd
from Bio import SeqIO
import gzip

def extract_metadata(fastq_file, output_csv):
    read_list = []
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
#        for record in SeqIO.parse(handle, "fastq"):
            desc_parts = record.description.split()
            meta = {p.split(":")[0]: p.split(":")[1] for p in desc_parts if ":" in p}
            if not all(k in meta for k in ["Start", "End", "Strand", "Chromosome", "CDS_info", "Source"]):
                continue
            try:
                start = int(meta["Start"])
                end = int(meta["End"])
            except:
                continue
            read_list.append({
                "read_idx": i,
                "read_id": record.id,
                "start": start,
                "end": end,
                "strand": meta["Strand"],
                "chromosome": meta["Chromosome"],
                "cds_info": meta["CDS_info"],
                "coding_status": meta["CDS_info"].split("_")[-1],
                "source": meta["Source"]
            })

    df = pd.DataFrame(read_list)
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} reads to {output_csv}")

if __name__ == "__main__":
    for i in range(1,5):
        extract_metadata(
            f"/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_{i}_128bact_128euk_R1.fastq",
            f"/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/mix_{i}_128bact_128euk_R1.csv"
        )
