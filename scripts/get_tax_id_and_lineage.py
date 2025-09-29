import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time

# Function to get Taxonomy ID from RefSeq Assembly Accession
def get_taxid_from_refseq(refseq_id):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "assembly",
        "term": refseq_id,
        "retmode": "json"
    }
    
    response = requests.get(base_url, params=params)
    data = response.json()
    
    if "esearchresult" in data and "idlist" in data["esearchresult"] and data["esearchresult"]["idlist"]:
        assembly_id = data["esearchresult"]["idlist"][0]  # First match
        return get_taxid_from_assembly(assembly_id)
    else:
        return None

# Function to get Taxonomy ID from Assembly ID
def get_taxid_from_assembly(assembly_id):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "assembly",
        "id": assembly_id,
        "retmode": "json"
    }

    response = requests.get(base_url, params=params)
    data = response.json()
    
    if "result" in data and assembly_id in data["result"]:
        return data["result"][assembly_id].get("taxid")
    else:
        return None

# Function to get lineage from Taxonomy ID using efetch
def get_lineage_from_taxid(taxid):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml"
    }

    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        return {}

    # Parse XML response
    root = ET.fromstring(response.content)
    lineage_ranks = {"phylum": None, "class": None, "order": None, "family": None}

    for taxon in root.findall(".//Taxon"):
        rank = taxon.find("Rank")
        name = taxon.find("ScientificName")

        if rank is not None and name is not None:
            rank_lower = rank.text.lower()
            if rank_lower in lineage_ranks:
                lineage_ranks[rank_lower] = name.text

    return lineage_ranks

input_csv = "/cs/usr/stavperez/sp/InSilicoSeq/InSilicoSeq/euk_org_info.csv"
df = pd.read_csv(input_csv)

# Add new columns for taxonomy data
df["Taxonomy ID"] = None
df["Phylum"] = None
df["Class"] = None
df["Order"] = None
df["Family"] = None

# Process each row
for index, row in df.iterrows():
    refseq_id = row["annotation"]
    taxid = get_taxid_from_refseq(refseq_id)

    if taxid:
        lineage = get_lineage_from_taxid(taxid)
    else:
        lineage = {}

    # Update DataFrame with results
    df.at[index, "Taxonomy ID"] = taxid
    df.at[index, "Phylum"] = lineage.get("phylum")
    df.at[index, "Class"] = lineage.get("class")
    df.at[index, "Order"] = lineage.get("order")
    df.at[index, "Family"] = lineage.get("family")

    print(f"Processed {refseq_id} -> Taxonomy ID: {taxid}, Lineage: {lineage}")
    time.sleep(1)  # To avoid API rate limits

# Save updated CSV
output_csv = "/cs/usr/stavperez/sp/InSilicoSeq/InSilicoSeq/euk_org_full_info.csv"
df.to_csv(output_csv, index=False)
print(f"\nUpdated CSV saved as: {output_csv}")
