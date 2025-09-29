import pandas as pd

# Load the annotated metadata CSV with the in_train flag
df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/chromosome_assemblies_with_train_flag.csv")

# Step 1: Identify genus names that contain at least one organism in the training set
excluded_genus = df[df["in_train"] == True]["genus"].unique()

# Step 2: Filter out those genera
df_filtered = df[~df["genus"].isin(excluded_genus)]

# Step 3: Keep only one organism per remaining genus (e.g., first occurrence)
df_unique_genus = df_filtered.drop_duplicates(subset=["genus"])

# Step 4: Save to a new CSV
df_unique_genus.to_csv("/cs/usr/stavperez/sp/genomes_new/chromosome_assemblies_unique_genus_non_train.csv", index=False)
