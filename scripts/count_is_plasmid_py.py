import pandas as pd

# === Load annotated GO file ===
file_path = "/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/final_generalized_GO_mix_1_with_plasmid.csv"
df = pd.read_csv(file_path)

# === Filter only relevant rows ===
df = df[df["confusion_type"].isin(["TP", "FN"])]

# === Contingency Table ===
count_table = pd.crosstab(df["confusion_type"], df["is_plasmid"])

print("=== Contingency Table (Counts) ===")
print(count_table)
print()

# === Totals ===
total_FN = df["confusion_type"].value_counts().get("FN", 0)
total_TP = df["confusion_type"].value_counts().get("TP", 0)

total_is_plasmid_true = df["is_plasmid"].sum()
total_is_plasmid_false = (~df["is_plasmid"]).sum()

# === Get intersection counts ===
FN_true = count_table.loc["FN", True] if True in count_table.columns else 0
TP_true = count_table.loc["TP", True] if True in count_table.columns else 0

# === Compute percentages ===
total_true_plasmid = FN_true + TP_true
perc_FN_true = (FN_true / total_FN) * 100 if total_true_plasmid > 0 else 0
perc_TP_true = (TP_true / total_TP) * 100 if total_true_plasmid > 0 else 0

# === Print totals ===
print("=== Totals ===")
print(f"Total FN: {total_FN}")
print(f"Total TP: {total_TP}")
print(f"Total is_plasmid == True: {total_is_plasmid_true}")
print(f"Total is_plasmid == False: {total_is_plasmid_false}")
print()

# === Print percentages with breakdown ===
print("=== Percentage Breakdown (only where is_plasmid == True) ===")
print(f"FN & is_plasmid == True: {FN_true} / {total_FN} = {perc_FN_true:.2f}%")
print(f"TP & is_plasmid == True: {TP_true} / {total_TP} = {perc_TP_true:.2f}%")
