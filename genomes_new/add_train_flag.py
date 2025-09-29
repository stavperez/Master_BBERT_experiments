import pandas as pd

# Load the bacteria metadata CSV
#meta_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/chromosome_assemblies_with_taxonomy.csv")
meta_df = pd.read_csv("/cs/usr/stavperez/sp/genomes_new/genome_refseq_bacteria.csv")

# Load the training abundance file and extract GCF/GCA IDs from the file names
train_df = pd.read_csv(
    "/sci/labs/aerez/aerez/backup/BBERTooD/training_data_only_bacteria/diverse_bacteria_training_abundance.txt",
    sep="\t", header=None, names=["fna_file", "abundance"]
)
train_df["assembly_accession"] = train_df["fna_file"].str.extract(r"^(GCF_\d+\.\d+|GCA_\d+\.\d+)")

# Add the "in_train" column by checking if assembly_accession is in the extracted list
meta_df["in_train"] = meta_df["assembly_accession"].isin(train_df["assembly_accession"])

# Save the updated CSV
meta_df.to_csv("/cs/usr/stavperez/sp/genomes_new/genome_refseq_bacteria_with_train_flag.csv", index=False)
