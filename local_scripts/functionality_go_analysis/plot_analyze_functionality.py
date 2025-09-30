import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu

# === Input and Output ===
input_dir = "volcano_plot_data_all"
output_dir = "volcano_plots_foldchange_fisher_mannwhitney"
os.makedirs(output_dir, exist_ok=True)

files = [
    "coding_status_volcano_data.csv",
    "GO_names_generalized_volcano_data.csv",
    "GO_slim_categories_volcano_data.csv",
    "GO_terms_volcano_data.csv",
    "product_volcano_data.csv"
]

# === Function to compute p-values and plot ===
def compute_stats_and_plot(df, label):
    df = df.copy()

    # === Fisher's Exact Test ===
    tp_total = df['TP'].sum()
    fn_total = df['FN'].sum()
    pvals_fisher = []
    for _, row in df.iterrows():
        fn = int(row['FN'])
        tp = int(row['TP'])
        rest_fn = fn_total - fn
        rest_tp = tp_total - tp
        table = [[fn, rest_fn], [tp, rest_tp]]
        try:
            _, p = fisher_exact(table)
        except:
            p = 1.0
        pvals_fisher.append(p)
    df['p_value_fisher'] = pvals_fisher
    df['-log10_pval_fisher'] = -np.log10(df['p_value_fisher'].replace(0, 1e-300))

    # === Mann–Whitney U Test ===
    pvals_mw = []
    for _, row in df.iterrows():
        fn = int(row['FN'])
        tp = int(row['TP'])
        if fn == 0 or tp == 0:
            pvals_mw.append(1.0)
        else:
            try:
                _, p = mannwhitneyu([1]*fn, [0]*tp, alternative="two-sided")
                pvals_mw.append(p)
            except:
                pvals_mw.append(1.0)
    df['p_value_mw'] = pvals_mw
    df['-log10_pval_mw'] = -np.log10(df['p_value_mw'].replace(0, 1e-300))

    # === Plot 1: Fold Change vs Fisher
    plt.figure(figsize=(10, 6))
    plt.scatter(df['log2_fold_change'], df['-log10_pval_fisher'], alpha=0.6, c='darkblue')
    plt.axhline(-np.log10(0.05), linestyle='--', color='black')
    plt.axvline(1, linestyle='--', color='gray')
    plt.axvline(-1, linestyle='--', color='gray')
    plt.title(f"Volcano Plot (Fisher) — {label}")
    plt.xlabel("log2(FN_norm / TP_norm)")
    plt.ylabel("-log10(p-value) [Fisher]")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{label}_volcano_fisher.png"), dpi=300)
    plt.close()

    # === Plot 2: Fold Change vs Mann–Whitney
    plt.figure(figsize=(10, 6))
    plt.scatter(df['log2_fold_change'], df['-log10_pval_mw'], alpha=0.6, c='darkred')
    plt.axhline(-np.log10(0.05), linestyle='--', color='black')
    plt.axvline(1, linestyle='--', color='gray')
    plt.axvline(-1, linestyle='--', color='gray')
    plt.title(f"Volcano Plot (Mann–Whitney) — {label}")
    plt.xlabel("log2(FN_norm / TP_norm)")
    plt.ylabel("-log10(p-value) [Mann–Whitney]")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{label}_volcano_mannwhitney.png"), dpi=300)
    plt.close()

# === Run on all files ===
for file in files:
    filepath = os.path.join(input_dir, file)
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
    df = pd.read_csv(filepath, index_col=0)
    label = file.replace("_volcano_data.csv", "")
    compute_stats_and_plot(df, label)

print(f"✅ Volcano plots saved to: {os.path.abspath(output_dir)}")
