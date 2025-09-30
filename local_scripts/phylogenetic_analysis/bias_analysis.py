import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, probplot, kruskal
import scikit_posthocs as sp


prefix = "unique_genus"
taxonomic_levels = ["phylum", "class", "order", "family"]


def check_normality(df, level):
    # Clean column names
    df.columns = [col.strip() for col in df.columns]
    df["FNR (%)"] = pd.to_numeric(df["FNR (%)"], errors="coerce")
    df = df.dropna(subset=["FNR (%)"])
    print(f"\n--- Normality test by {level.upper()} ---")
    grouped = df.groupby(level)

    for group_name, group_df in grouped:
        values = group_df["FNR (%)"].dropna()

        if len(values) < 3:
            continue  # not enough data

        stat, p = shapiro(values)
        status = "NOT normal" if p < 0.05 else "Normal"
        print(
            f"{group_name:<30}  N={len(values):<3}  p={p:.4f}  {status}")
    sns.boxplot(x=level, y="FNR (%)", data=df)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def Kruskal_Wallis(df, level):
    df.columns = [col.strip() for col in df.columns]
    df["FNR (%)"] = pd.to_numeric(df["FNR (%)"], errors="coerce")
    df = df.dropna(subset=["FNR (%)"])
    print(f"\n--- Kruskal–Wallis Test by {level.upper()} ---")

    # Group values by taxonomic level
    grouped_values = [
        group["FNR (%)"].dropna().values
        for name, group in df.groupby(level)
        if len(group) >= 3  # only keep groups with 3+ samples
    ]

    group_names = [
        name for name, group in df.groupby(level)
        if len(group) >= 3
    ]

    stat, p = kruskal(*grouped_values)
    print(f"Kruskal–Wallis H = {stat:.4f}, p = {p:.4g}")

    if p < 0.05:
        print("=> Significant differences found between groups.")
    else:
        print("=> No significant differences found between groups.")

    # Boxplot
    plt.figure(figsize=(12, 5))
    sns.boxplot(data=df[df[level].isin(group_names)],
                x=level, y="FNR (%)")
    plt.title(f"FNR (%) by {level}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


import scikit_posthocs as sp

def run_dunn_and_find_biased_groups(df, level):
    print(f"\nRunning Dunn’s test and detecting high-bias groups at level: {level.upper()}")

    # Filter groups with enough samples
    filtered_df = df.groupby(level).filter(lambda x: len(x) >= 3)
    if filtered_df[level].nunique() < 2:
        print("Not enough valid groups with ≥ 3 samples.")
        return

    # Dunn’s test
    dunn_results = sp.posthoc_dunn(
        filtered_df, val_col="FNR (%)", group_col=level, p_adjust="holm"
    )

    # Compute mean FNR per group
    group_means = filtered_df.groupby(level)["FNR (%)"].mean().sort_values(ascending=False)

    # Detect biased groups
    biased_groups = []

    for group in group_means.index:
        significant_vs = []
        for other in group_means.index:
            if group == other:
                continue
            pval = dunn_results.loc[group, other] if pd.notna(dunn_results.loc[group, other]) else dunn_results.loc[other, group]
            if pval < 0.05 and group_means[group] > group_means[other]:
                significant_vs.append((other, pval))

        if len(significant_vs) >= 2:
            biased_groups.append((group, group_means[group], significant_vs))

    # Print result
    if biased_groups:
        print("\nBiased groups with significantly higher FNR than ≥2 others:")
        for g, mean, comparisons in biased_groups:
            print(f"- {g} (mean FNR = {mean:.2f}%) > {len(comparisons)} groups:")
            for other, pval in sorted(comparisons, key=lambda x: x[1]):
                print(f"   - vs {other}: p = {pval:.4g}")
    else:
        print("No strongly biased groups found.")

def plot_group_distributions(df, level):
    filtered_df = df.groupby(level).filter(lambda x: len(x) >= 3)
    if filtered_df.empty:
        print("No groups with enough data to plot.")
        return

    plt.figure(figsize=(12, 5))
    sns.boxplot(data=filtered_df, x=level, y="FNR (%)")
    plt.title(f"FNR (%) by {level}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def main():
    df = pd.read_csv("FNR_all.csv")
    for level in taxonomic_levels:
        # Check normality
        check_normality(df, level)
        # Check if there are differences
        Kruskal_Wallis(df, level)
        # Check differences between each 2 groups
        run_dunn_and_find_biased_groups(df, level)
        # Plot
        plot_group_distributions(df, level)


if __name__ == "__main__":
    main()
