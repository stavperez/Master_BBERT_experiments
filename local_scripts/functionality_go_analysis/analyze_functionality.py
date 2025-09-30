import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

print("--- Starting Data Analysis Script ---")

# === Load the data ===
input_file = r"C:\Users\stavp\CSE\Master\Lab\Project\plot_embeddings\data" \
             r"\Pseudomonas_Saccharomyces_GO_info.csv"

output_dir = "volcano_plot_data_Pseudomonas_Saccharomyces"
os.makedirs(output_dir, exist_ok=True)

functional_columns = [
    "GO_names_generalized",
    "product",
    "GO_terms",
    "coding_status",
    "GO_slim_categories"
]

# === Load CSV once ===
print(f"\n1. Loading input CSV: {input_file}")
df = pd.read_csv(input_file)
print(f"   - Total rows: {len(df)}")

# === Filter for relevant rows ===
df = df[df['product'].notna() | df['GO_terms'].notna()]
print(f"   - After filtering: {len(df)} rows with 'product' or 'GO_terms'")

# === Only keep columns we need ===
columns_to_keep = [
    'read_id', 'confusion_type', 'cds_info', 'gene', 'functional_annotation',
    'product', 'GO_terms', 'GO_terms_generalized', 'GO_names_generalized',
    'GO_slim_categories', 'cds_overlap_pct', 'coding_status'
]
df = df[columns_to_keep]

# === Function to process a single functional column ===
def process_column(df, column_name):
    print(f"\n--- Processing column: {column_name} ---")

    # Explode semicolon-delimited values
    df_expanded = df[df['confusion_type'].isin(['TP', 'FN'])].copy()
    df_expanded[column_name] = df_expanded[column_name].fillna('').astype(str)
    df_expanded[column_name] = df_expanded[column_name].str.split(';')
    df_expanded = df_expanded.explode(column_name)
    df_expanded[column_name] = df_expanded[column_name].str.strip()
    df_expanded = df_expanded[df_expanded[column_name] != '']

    # Count TP/FN per term
    print("   - Counting TP/FN per term...")
    term_stats = df_expanded.groupby([column_name, 'confusion_type']).size().unstack(fill_value=0)
    term_stats['TP'] = term_stats.get('TP', 0)
    term_stats['FN'] = term_stats.get('FN', 0)

    # Global totals
    tp_total = term_stats['TP'].sum()
    fn_total = term_stats['FN'].sum()

    print(f"   - Total TP reads: {tp_total}")
    print(f"   - Total FN reads: {fn_total}")

    # Normalized proportions
    term_stats['TP_norm'] = term_stats['TP'] / tp_total
    term_stats['FN_norm'] = term_stats['FN'] / fn_total

    # Log2 fold change (with pseudocount)
    term_stats['log2_fold_change'] = np.log2((term_stats['FN_norm'] + 1e-10) / (term_stats['TP_norm'] + 1e-10))

    # Fisher exact test for significance
    print("   - Running Fisher’s exact test for each term...")
    pvals = []
    for idx, row in term_stats.iterrows():
        fn = row['FN']
        tp = row['TP']
        fn_other = fn_total - fn
        tp_other = tp_total - tp
        contingency = [[fn, fn_other], [tp, tp_other]]
        try:
            _, p = fisher_exact(contingency)
        except:
            p = 1.0
        pvals.append(p)

    term_stats['p_value'] = pvals
    term_stats['-log10_p_value'] = -np.log10(term_stats['p_value'].replace(0, 1e-300))

    # Save table
    result_df = term_stats[['FN', 'TP', 'FN_norm', 'TP_norm', 'log2_fold_change', 'p_value', '-log10_p_value']]
    result_df.to_csv(os.path.join(output_dir, f"{column_name}_volcano_data.csv"))
    print(f"   - Saved results to: {column_name}_volcano_data.csv")

    return result_df

# === Volcano plot function ===
def plot_volcano(df, column_name):
    print(f"   - Plotting volcano plot for: {column_name}")
    df = df.copy()
    df['category'] = 'Neither'
    df.loc[(df['log2_fold_change'] > 1) & (df['p_value'] < 0.05), 'category'] = 'FN-Enriched'
    df.loc[(df['log2_fold_change'] < -1) & (df['p_value'] < 0.05), 'category'] = 'TP-Enriched'

    plt.figure(figsize=(10, 7))
    palette = {'FN-Enriched': 'red', 'TP-Enriched': 'blue', 'Neither': 'gray'}
    for cat, color in palette.items():
        subset = df[df['category'] == cat]
        plt.scatter(subset['log2_fold_change'], subset['-log10_p_value'], c=color, label=cat, s=10, alpha=0.6)

    plt.axhline(-np.log10(0.05), linestyle='--', color='black')
    plt.axvline(1, linestyle='--', color='black')
    plt.axvline(-1, linestyle='--', color='black')
    plt.xlabel('log2(FN_norm / TP_norm)')
    plt.ylabel('-log10(p-value)')
    plt.title(f'Volcano Plot: {column_name}')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)
    plot_path = os.path.join(output_dir, f"{column_name}_volcano_plot.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"   - Saved plot to: {plot_path}")

# === Process all functional columns ===
for col in functional_columns:
    result = process_column(df, col)
    plot_volcano(result, col)

print("\n--- DONE ---")
