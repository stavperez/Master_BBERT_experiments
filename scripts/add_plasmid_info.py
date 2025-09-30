import pandas as pd

# === Input files ===
go_file = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/final_generalized_GO_mix_1.csv"
plasmid_map_file = "/cs/usr/stavperez/sp/genomes_new/sequence_header_to_plasmid.csv"
output_file = go_file.replace(".csv", "_with_plasmid.csv")

# === Load data ===
df_go = pd.read_csv(go_file)
df_plasmid = pd.read_csv(plasmid_map_file)

# === Build nested dict: {organism_name: {sequence_header: is_plasmid}} ===
plasmid_lookup = {}
for _, row in df_plasmid.iterrows():
    org = row["organism_name"]
    header = row["sequence_header"]
    is_plasmid = row["is_plasmid"]

    if org not in plasmid_lookup:
        plasmid_lookup[org] = {}
    plasmid_lookup[org][header] = is_plasmid

# === Annotate each row in GO file ===
def get_is_plasmid(row):
    org = row["source_x"]
    chrom = row["chromosome"]
    return plasmid_lookup.get(org, {}).get(chrom, False)

df_go["is_plasmid"] = df_go.apply(get_is_plasmid, axis=1)

# === Save result ===
df_go.to_csv(output_file, index=False)
print(f"[DONE] Saved annotated file with 'is_plasmid' column to: {output_file}")
