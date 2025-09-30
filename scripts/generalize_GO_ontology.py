import pandas as pd
from goatools.obo_parser import GODag

# === CONFIGURATION ===
INPUT_CSV = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/final_annotated_reads_mix_1.csv"
OUTPUT_CSV = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/final_generalized_GO_mix_1.csv"
#INPUT_CSV = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/annotated_reads_split_mix_1/Pseudomonas_aeruginosa_annotated.csv"
#OUTPUT_CSV = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/embeddings/output_scores/Pseudomonas_aeruginosa/Pseudomonas_aeruginosa_go_terms.csv"
FULL_OBO = "go-basic.obo"
GOSLIM_OBO = "goslim_generic.obo"
MAX_DEPTH = 3

# === Load Ontologies ===
print("Loading GO full ontology...")
go_dag = GODag(FULL_OBO)

print("Loading GO Slim (generic)...")
go_slim = GODag(GOSLIM_OBO)
go_slim_ids = set(go_slim.keys())

# === Helper Functions ===
def get_top_level_terms(go_id, max_depth=2):
    if go_id not in go_dag:
        return set()
    term = go_dag[go_id]
    result = set()
    def climb(node, depth):
        if depth > max_depth:
            return
        result.add(node.id)
        for parent in node.parents:
            climb(parent, depth + 1)
    climb(term, 0)
    return result

def parse_and_generalize(go_str):
    if pd.isna(go_str) or go_str == "Non-CDS":
        return ""
    ids = [p.split("|")[1] for p in go_str.split(",") if "|" in p]
    ids = [f"GO:{p.zfill(7)}" for p in ids if p.isdigit()]
    generalized = set()
    for gid in ids:
        if gid in go_dag:
            generalized.update(get_top_level_terms(gid, MAX_DEPTH))
    return ";".join(sorted(generalized))

def go_to_name(go_id):
    return go_dag[go_id].name if go_id in go_dag else "UNKNOWN"

def map_to_slim(go_id):
    if go_id not in go_dag:
        return None
    for ancestor_id in go_dag[go_id].get_all_parents():
        if ancestor_id in go_slim_ids:
            return go_dag[ancestor_id].name  
    return None

def names_from_ids(ids_str):
    return ";".join(go_to_name(gid) for gid in ids_str.split(";") if gid.strip())

def slim_categories(ids_str):
    cats = {map_to_slim(gid) for gid in ids_str.split(";") if gid.strip()}
    return ";".join(sorted(cat for cat in cats if cat))

# === Process ===
print("Reading input CSV...")
df = pd.read_csv(INPUT_CSV)
print(f"length of input df: {len(df)} reads")

print("Generalizing GO terms using the full hierarchy...")
df["GO_terms_generalized"] = df["GO_terms"].apply(parse_and_generalize)

print("Adding human-readable GO names...")
df["GO_names_generalized"] = df["GO_terms_generalized"].apply(names_from_ids)

print("Mapping to GO Slim categories...")
df["GO_slim_categories"] = df["GO_terms_generalized"].apply(slim_categories)

print("Saving enhanced CSV...")
df.to_csv(OUTPUT_CSV, index=False)
print(OUTPUT_CSV)
print("Done.")
