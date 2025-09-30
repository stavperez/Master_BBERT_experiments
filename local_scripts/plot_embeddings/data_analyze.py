import pandas as pd

# Load data
df = pd.read_csv("data/Pseudomonas_Saccharomyces_GO_info_with_frame.csv")

# Normalize reading_frame to include None as string
df["reading_frame"] = df["reading_frame"].fillna("None").astype(str)

# Create 3D table (CDS x confusion_type x reading_frame)
pivot_table = pd.crosstab(
    index=[df["CDS"], df["confusion_type"]],
    columns=df["reading_frame"],
    dropna=False
)

print("\n=== 3D Table (CDS x confusion_type x reading_frame) ===")
print(pivot_table)

# Calculate FN / (FN + TP) for coding & None frame
fn_count = len(df[(df["CDS"] == "coding") & (df["reading_frame"] == "None") & (df["confusion_type"] == "FN")])
tp_count = len(df[(df["CDS"] == "coding") & (df["reading_frame"] == "None") & (df["confusion_type"] == "TP")])
total = fn_count + tp_count

fn_total = len(df[(df["CDS"] == "coding") & (df["confusion_type"] == "FN")])
tp_total = len(df[(df["CDS"] == "coding") & (df["confusion_type"] == "TP")])


print("\n=== FN vs TP for coding & reading_frame=None ===")
print(f"FN: {fn_count}, TP: {tp_count}, Total FN: {fn_total}, Total TP: {tp_total}")
print(f"Percentage FN: {(100 * fn_count / fn_total):.2f}%")
print(f"Percentage TP: {(100 * tp_count / tp_total):.2f}%")

