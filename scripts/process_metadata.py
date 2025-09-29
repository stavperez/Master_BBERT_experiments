import json
import pandas as pd
from ete3 import NCBITaxa

# Optionally, set the custom cache location if not done via environment variable:
# config.ncbi_taxa_cache = "/sci/labs/aerez/stavperez/.etetoolkit"

ncbi = NCBITaxa()

def get_lineage_info(tax_id):
    try:
        lineage = ncbi.get_lineage(tax_id)
        names = ncbi.get_taxid_translator(lineage)
        # Create a string representation of the lineage
        lineage_str = "; ".join([names[t] for t in lineage])
        return lineage_str
    except Exception as e:
        print(f"Error retrieving lineage for tax_id {tax_id}: {e}")
        return "Unknown"

#jsonl_file = "eukaryotes_metadata_lineage.jsonl"
#jsonl_file = "prokaryotes_metadata_lineage.jsonl"
organism_data = []

with open(jsonl_file, "r") as f:
    for line in f:
        record = json.loads(line)
        organism_name = record.get("organism", {}).get("organism_name", "Unknown")
        tax_id = record.get("organism", {}).get("tax_id", "Unknown")
        refseq_category = record.get("assembly", {}).get("refseq_category", "Unknown")
        
        # Only try to get lineage if tax_id is valid
        if tax_id != "Unknown":
            lineage = get_lineage_info(tax_id)
        else:
            lineage = "Unknown"

        organism_data.append({
            "Organism Name": organism_name,
            "Tax ID": tax_id,
            "RefSeq Category": refseq_category,
            "Lineage": lineage
        })

df = pd.DataFrame(organism_data)
print("Number of records in DataFrame:", len(df))
df.to_csv("prokaryotes_organisms_with_genomes.csv", index=False)
print("Saved prokaryotic organisms with RefSeq genomes and lineage information to 'prokaryotic_organisms_with_genomes.csv'")
