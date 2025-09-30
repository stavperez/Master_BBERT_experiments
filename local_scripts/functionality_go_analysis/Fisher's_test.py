import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import glob
import os

# ----------- CONFIGURATION -------------
TOP_LABELS = 10  # Top N terms to label in volcano plot
TOP_BARS = 20    # Top N enriched terms to show in bar plot
PVAL_THRESHOLD = 0.05
INPUT_PATTERN = "*volcano_data.csv"
# ---------------------------------------

def process_file(file_path):
    print(f"\nProcessing: {file_path}")
    df = pd.read_csv(file_path)

    # Infer the name of the "term" column from the file name
    filename = os.path.basename(file_path)
    term_col_candidate = filename.replace('_volcano_data.csv', '')
    if term_col_candidate not in df.columns:
        raise ValueError(f"Expected column '{term_col_candidate}' not found in {file_path}")

    # Rename it to 'term' for consistency
    df = df.rename(columns={term_col_candidate: 'term'})

    # Check required columns
    required_cols = {'term', 'TP', 'FN', 'TP_norm', 'FN_norm'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in {file_path}: {required_cols - set(df.columns)}")

    # Compute totals
    total_TP = df['TP'].sum()
    total_FN = df['FN'].sum()

    fisher_pvals = []
    odds_ratios = []

    for _, row in df.iterrows():
        TP = int(row['TP'])
        FN = int(row['FN'])
        not_TP = total_TP - TP
        not_FN = total_FN - FN
        table = [[TP, FN], [not_TP, not_FN]]
        odds_ratio, p = fisher_exact(table, alternative='two-sided')
        fisher_pvals.append(p)
        odds_ratios.append(odds_ratio)

    # Add results
    df['fisher_p'] = fisher_pvals
    df['odds_ratio'] = odds_ratios

    # If log2_fold_change is missing, compute it
    if 'log2_fold_change' not in df.columns:
        df['log2_fold_change'] = np.log2((df['TP_norm'] + 1e-9) / (df['FN_norm'] + 1e-9))

    # Save updated file
    out_file = file_path.replace('.csv', '_with_fisher.csv')
    df.to_csv(out_file, index=False)
    print(f"Saved enriched data to: {out_file}")

    # ---- Volcano Plot ----
    df['neglog10_p'] = -np.log10(df['fisher_p'] + 1e-300)

    plt.figure(figsize=(10, 6))
    plt.scatter(
        df['log2_fold_change'],
        df['neglog10_p'],
        s=10,  # point size (default is ~36)
        alpha=0.6,  # transparency (0=fully transparent, 1=opaque)
        color='gray',
        label='All terms'
    )

    # Significant terms: slightly larger and red
    sig = df['fisher_p'] < PVAL_THRESHOLD
    plt.scatter(
        df[sig]['log2_fold_change'],
        df[sig]['neglog10_p'],
        s=20,
        alpha=0.3,
        color='red',
        label='p < 0.05'
    )
    # top_sig = df[sig].nsmallest(TOP_LABELS, 'fisher_p')
    # for _, row in top_sig.iterrows():
    #     plt.text(row['log2_fold_change'], row['neglog10_p'], row['term'], fontsize=8)

    # Axis lines and labels
    plt.axhline(-np.log10(PVAL_THRESHOLD), color='gray', linestyle='--',
                linewidth=1)
    plt.axvline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('log₂ fold change (TP_norm / FN_norm)')
    plt.ylabel('-log₁₀(p-value)')
    plt.title('Volcano Plot')
    plt.legend()
    plt.tight_layout()

    volcano_plot_path = file_path.replace('.csv', '_volcano_plot.png')
    plt.savefig(volcano_plot_path)
    print(f"Saved volcano plot to: {volcano_plot_path}")
    plt.show()
    # ---- Bar Plot of TP_norm vs FN_norm ----
    top_enriched = df[sig].sort_values('log2_fold_change', ascending=False).head(TOP_BARS)
    bar_df = top_enriched[['term', 'TP_norm', 'FN_norm']].set_index('term')

    bar_df.plot(kind='bar', figsize=(12, 6), title='Normalized Frequency: TP_norm vs FN_norm (Top Enriched Terms)', rot=45)
    plt.ylabel('Normalized frequency')
    plt.tight_layout()

    bar_plot_path = file_path.replace('.csv', '_bar_plot.png')
    plt.savefig(bar_plot_path)
    print(f"Saved bar plot to: {bar_plot_path}")
    plt.show()

    ### histogram

    # Optional: Histogram of odds_ratio (significant terms only)
    sig_df = df[df['fisher_p'] < PVAL_THRESHOLD].copy()
    sig_df = sig_df.replace([np.inf, 0],
                            np.nan)  # Remove inf and 0 to avoid issues

    # Drop missing values
    sig_df = sig_df.dropna(subset=['odds_ratio'])

    # Compute log2(odds_ratio)
    sig_df['log2_odds_ratio'] = np.log2(sig_df['odds_ratio'])

    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.hist(sig_df['log2_odds_ratio'], bins=30, color='skyblue',
             edgecolor='black')
    plt.axvline(0, color='black', linestyle='--')
    plt.title(
        'Distribution of log₂(odds ratio) for significant terms (p < 0.05)')
    plt.xlabel('log₂(odds ratio)')
    plt.ylabel('Number of terms')
    plt.tight_layout()

    hist_path = file_path.replace('.csv', '_log2_odds_histogram.png')
    plt.savefig(hist_path)
    print(f"Saved odds ratio histogram to: {hist_path}")
    plt.show()


# ---------- Run on all matching files ----------
for csv_file in glob.glob(INPUT_PATTERN):
    process_file(csv_file)
