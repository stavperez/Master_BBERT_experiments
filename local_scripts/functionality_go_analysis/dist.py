import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# === List of files and their associated go_term_category column ===
files_and_columns = {
    "coding_status_volcano_data_with_fisher.csv": "coding_status",
    "GO_names_generalized_volcano_data_with_fisher.csv": "GO_names_generalized",
    "GO_slim_categories_volcano_data_with_fisher.csv": "GO_slim_categories",
    "GO_terms_volcano_data_with_fisher.csv": "GO_terms",
    "product_volcano_data_with_fisher.csv": "product"
}

# === Output directory for plots ===
output_dir = "fisher_plots"
os.makedirs(output_dir, exist_ok=True)

# === Process each file ===
for file_path, term_column in files_and_columns.items():
    print(f"\n=== Processing: {file_path} ===")

    # Load CSV
    df = pd.read_csv(file_path)

    if term_column not in df.columns:
        raise ValueError(f"Column '{term_column}' not found in {file_path}")

    # Recompute log2_fold_change with pseudo-count 0.5
    df["log2_fold_change"] = np.log2((df["FN_norm"] + 0.5) / (df["TP_norm"] + 0.5))

    # Plot distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df["log2_fold_change"], bins=50, edgecolor="black")
    plt.title(f"log2_fold_change Distribution: {term_column}")
    plt.xlabel("log2_fold_change")
    plt.ylabel("Count")
    plt.grid(True)
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f"{term_column}_distribution.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    # Compute quantiles
    lower_thresh = df["log2_fold_change"].quantile(0.025)
    upper_thresh = df["log2_fold_change"].quantile(0.975)

    # Filter top and bottom
    low_terms = df[df["log2_fold_change"] <= lower_thresh]
    high_terms = df[df["log2_fold_change"] >= upper_thresh]

    # Print results with actual category values
    print(f"\n--- Bottom 2.5% log2_fold_change for {term_column} ---")
    print(low_terms[[term_column, "log2_fold_change"]].sort_values("log2_fold_change"))

    print(f"\n--- Top 2.5% log2_fold_change for {term_column} ---")
    print(high_terms[[term_column, "log2_fold_change"]].sort_values("log2_fold_change", ascending=False))
