import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score



prefix = "unique_genus"
# prefix = "new"

def scatter_plot(category, df, ds=1, size=128):
    global min_val, max_val
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x=f"Full_{category} (%)",
                    y=f"Classified_{category} (%)")
                    # , label="Genomes")
    # Add x = y line
    min_val = min(df[f"Full_{category} (%)"].min(),
                  df[f"Classified_{category} (%)"].min())
    max_val = max(df[f"Full_{category} (%)"].max(),
                  df[f"Classified_{category} (%)"].max())
    plt.plot([min_val, max_val], [min_val, max_val], color='red',
             linestyle='--', linewidth=1.5, label='x = y')
    # Fit linear regression
    x = df[f"Full_{category} (%)"].values.reshape(-1, 1)
    y = df[f"Classified_{category} (%)"].values
    reg = LinearRegression(fit_intercept=False).fit(x, y)
    y_pred = reg.predict(x)
    # Plot regression line
    plt.plot(df[f"Full_{category} (%)"], y_pred, color='blue', linewidth=2,
             label="Linear Fit")
    # Calculate regression parameters
    slope = reg.coef_[0]
    intercept = reg.intercept_
    r2 = r2_score(y, y_pred)
    # Add regression equation and R² to plot
    eq_text = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r2:.3f}"
    plt.text(0.05, 0.95, eq_text, transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top',
             bbox=dict(boxstyle="round", alpha=0.2))
    plt.title(f"Full vs. Classified Genome {category}\nSample size {size} "
              f"dataset "
              f"{ds}", fontsize=22)
    plt.xlabel(f"Full {category} (%)", fontsize=18)
    plt.ylabel(f"Classified {category} (%)", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend().remove()
    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_{category}_{ds}_{size}.png")
    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_{category}_"
                rf"{ds}_{size}.svg")
    plt.show()


def scatter_plot1(category, df, ds=1, size=128):
    global min_val, max_val

    # Determine correct columns based on category
    if category == "Coverage":
        x_col_name = f"Full_{category} (%)"
        y_col_name = f"Classified_{category} (%)"
    else:
        x_col_name = f"Full_{category}_from_covered (%)"
        y_col_name = f"Classified_{category}_from_covered (%)"

    plt.figure(figsize=(8, 6))
    # sns.scatterplot(data=df, x=x_col_name, y=y_col_name)
    # Plot scatter and regression with confidence interval
    sns.regplot(
        data=df, x=x_col_name, y=y_col_name,
        ci=95, scatter_kws={'s': 50}, line_kws={'color': 'blue'}, truncate=False
    )
    # Add x = y line
    min_val = min(df[x_col_name].min(), df[y_col_name].min())
    max_val = max(df[x_col_name].max(), df[y_col_name].max())
    plt.plot([min_val, max_val], [min_val, max_val], color='red',
             linestyle='--', linewidth=1.5, label='x = y')

    # Fit linear regression
    x = df[x_col_name].values.reshape(-1, 1)
    y = df[y_col_name].values
    reg = LinearRegression(fit_intercept=False).fit(x, y)
    y_pred = reg.predict(x)

    # Plot regression line
    plt.plot(df[x_col_name], y_pred, color='blue', linewidth=2, label="Linear Fit")

    # Calculate regression parameters
    slope = reg.coef_[0]
    intercept = reg.intercept_
    r2 = r2_score(y, y_pred)

    # Add regression equation and R² to plot
    eq_text = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r2:.3f}"
    plt.text(0.05, 0.95, eq_text, transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top',
             bbox=dict(boxstyle="round", alpha=0.2))

    plt.title(f"Full vs. Classified Genome {category}\nSample size {size} dataset {ds}",
              fontsize=22)
    plt.xlabel(f"Full {category} (%)", fontsize=18)
    plt.ylabel(f"Classified {category} (%)", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend().remove()

    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_{category}_{ds}_{size}.png")
    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_{category}_{ds}_{size}.svg")
    plt.show()

def scatter_plot2(category, df, ds=1, size=128):
    global min_val, max_val

    # Determine correct columns based on category
    if category == "Coverage":
        x_col_name = f"Full_{category} (%)"
        y_col_name = f"Classified_{category} (%)"
    else:
        x_col_name = f"Full_{category}_from_covered (%)"
        y_col_name = f"Classified_{category}_from_covered (%)"

    plt.figure(figsize=(8, 6))

    # Plot scatter and regression with 95% confidence interval
    sns.regplot(
        data=df, x=x_col_name, y=y_col_name,
        ci=95, scatter_kws={'s': 50}, line_kws={'color': 'blue'}, truncate=False
    )

    # Add x = y line
    min_val = min(df[x_col_name].min(), df[y_col_name].min())
    max_val = max(df[x_col_name].max(), df[y_col_name].max())
    plt.plot([min_val, max_val], [min_val, max_val], color='red',
             linestyle='--', linewidth=1.5, label='x = y')

    # Compute regression stats (for equation only)
    x = df[x_col_name].values.reshape(-1, 1)
    y = df[y_col_name].values
    reg = LinearRegression(fit_intercept=False).fit(x, y)
    slope = reg.coef_[0]
    intercept = reg.intercept_
    r2 = r2_score(y, reg.predict(x))

    # Add regression equation and R² to plot
    eq_text = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r2:.3f}"
    plt.text(0.05, 0.95, eq_text, transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top',
             bbox=dict(boxstyle="round", alpha=0.2))

    plt.title(f"Full vs. Classified Genome {category}\nSample size {size} dataset {ds}",
              fontsize=22)
    plt.xlabel(f"Full {category} (%)", fontsize=18)
    plt.ylabel(f"Classified {category} (%)", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend().remove()

    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_"
                rf"{category}_{ds}_{size}1.png")
    plt.savefig(rf"{prefix}\new_plots\scatter_full_vs_classified_"
                rf"{category}_{ds}_{size}1.svg")
    plt.show()

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def scatter_plot3(category, df):
    if category == "Coverage":
        x_col_name = f"Full_{category} (%)"
        y_col_name = f"Classified_{category} (%)"
    else:
        x_col_name = f"Full_{category}_from_covered (%)"
        y_col_name = f"Classified_{category}_from_covered (%)"

    # Map dataset codes to readable labels
    rename_map = {
        "1_128": "Dataset 1",
        "2_128": "Dataset 2",
        "3_128": "Dataset 3",
        "4_128": "Dataset 4"
    }
    df["ds_mapped"] = df["ds"].map(rename_map).fillna(df["ds"])

    # Create wider plot to allow space on the right
    fig, ax = plt.subplots(figsize=(10, 6))
    plt.subplots_adjust(right=0.7, top=0.9)

    # Scatter plot
    sns.scatterplot(
        data=df,
        x=x_col_name,
        y=y_col_name,
        hue="ds_mapped",
        palette="pastel",
        s=50,
        ax=ax,
        legend="full"
    )

    # Regression line
    sns.regplot(
        data=df,
        x=x_col_name,
        y=y_col_name,
        ci=95,
        scatter=False,
        line_kws={'color': 'blue'},
        truncate=False,
        ax=ax
    )
    ax.lines[-1].set_label("Regression line")

    # x = y reference line
    min_val = min(df[x_col_name].min(), df[y_col_name].min())
    max_val = max(df[x_col_name].max(), df[y_col_name].max())
    ax.plot([min_val, max_val], [min_val, max_val], color='red',
            linestyle='--', linewidth=2, label='x = y')

    # Compute regression stats
    x = df[x_col_name].values.reshape(-1, 1)
    y = df[y_col_name].values
    reg = LinearRegression(fit_intercept=False).fit(x, y)
    slope = reg.coef_[0]
    intercept = reg.intercept_
    r2 = r2_score(y, reg.predict(x))
    eq_text = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r2:.3f}"

    # Place equation outside plot on the right (above legend)
    fig.text(0.72, 0.88, eq_text, ha='left', va='top',
             fontsize=16, bbox=dict(boxstyle="round", alpha=0.2))

    # Axes formatting
    ax.set_title(f"Full vs. Classified Genome {category}", fontsize=22)
    ax.set_xlabel(f"Full {category} (%)", fontsize=18)
    ax.set_ylabel(f"Classified {category} (%)", fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True)

    # Move legend to the right, below equation
    legend = ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left',
                       title="Legend", fontsize=14, title_fontsize=16)
    for text in legend.get_texts():
        text.set_fontsize(14)

    # Save
    plt.savefig(rf"{prefix}/new_plots/scatter_full_vs_classified_{category}_all.png", dpi=300, bbox_inches='tight')
    plt.savefig(rf"{prefix}/new_plots/scatter_full_vs_classified_{category}_all.svg", bbox_inches='tight')
    plt.show()



def box_plot(category, df, ds=1, size=128):
    melted_df = df[[f"Full_{category} (%)", f"Classified_{category} (%)"]].melt(
        var_name=f"{category} Type", value_name="Percentage"
    )
    # Add x = y line
    min_val = min(df[f"Full_{category} (%)"].min(),
                  df[f"Classified_{category} (%)"].min())
    max_val = max(df[f"Full_{category} (%)"].max(),
                  df[f"Classified_{category} (%)"].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='x = y')

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=melted_df, x=f"{category} Type", y="Percentage")
    plt.title(
        f"{category} Distribution: Full vs Classified size {size} dataset {ds}")
    plt.xlabel(f"{category} Type")
    plt.ylabel(f"{category} (%)")
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.savefig(rf"{prefix}\plots\boxplot_full_vs_classified_{category}_{ds}_{size}.png")
    plt.show()

def plot_error_type_comparison(df, ds, size):
    error_types = ["%Single_Errors", "%Adjacent_Pairs", "%Longer_Runs"]
    plt.figure(figsize=(8, 6))

    colors = {"%Single_Errors": "red", "%Adjacent_Pairs": "orange", "%Longer_Runs": "green"}

    for et in error_types:
        full_col = f"Full_{et}"
        class_col = f"Classified_{et}"
        if full_col in df.columns and class_col in df.columns:
            plt.scatter(df[full_col], df[class_col], label=et.replace("_", " ").replace("%", ""), alpha=0.7, color=colors[et])

    min_val = 0
    max_val = max(df[[f"Full_{et}" for et in error_types] + [f"Classified_{et}" for et in error_types]].max())

    plt.plot([min_val, max_val], [min_val, max_val], 'k--', label="x = y")

    plt.title(f"Full vs. Classified Error Percentages\nSample size {size} "
              f"dataset "
              f"{ds}", fontsize=22)
    plt.xlabel("Full (%) Error Type", fontsize=18)
    plt.ylabel("Classified (%) Error Type", fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(rf"{prefix}\plots\scatter_error_type_percentages_{ds}_{size}.png")
    plt.savefig(rf"{prefix}\plots\scatter_error_type_percentages_{ds}_"
                rf"{size}.svg")
    # plt.show()

def plot_error_metrics(df, ds, size, mode):
    """
    mode: "of_Total" or "of_Genome"
    """
    error_types = ["Single_Errors", "Adjacent_Pairs", "Longer_Runs"]
    colors = {"Single_Errors": "red", "Adjacent_Pairs": "orange", "Longer_Runs": "green"}
    plt.figure(figsize=(8, 6))

    for et in error_types:
        full_col = f"Full_%{et}_{mode}"
        class_col = f"Classified_%{et}_{mode}"
        if full_col in df.columns and class_col in df.columns:
            plt.scatter(df[full_col], df[class_col], label=et.replace("_", " "), alpha=0.7, color=colors[et])

    all_vals = []
    for et in error_types:
        all_vals.extend(df[[f"Full_%{et}_{mode}", f"Classified_%{et}_{mode}"]].values.flatten())
    min_val = 0
    max_val = max(all_vals)

    plt.plot([min_val, max_val], [min_val, max_val], 'k--', label="x = y")
    plt.title(f"Error Type % Comparison ({mode})\nSample size {size} "
              f"dataset "
              f"{ds}", fontsize=22)
    plt.xlabel("Full (%)", fontsize=18)
    plt.ylabel("Classified (%)", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    mode_tag = mode.replace("of_", "")
    plt.savefig(rf"{prefix}\new_plots\scatter_error_type_{mode_tag}_{ds}_{size}.png")
    plt.savefig(rf"{prefix}\new_plots\scatter_error_type_{mode_tag}_{ds}_"
                rf"{size}.svg")
    # plt.show()

def plot_errors_distrubution(df, ds):
    # Plot the distribution
    sns.histplot(df['errors_outof_coverage_full_dist (%)'], bins=100,
                 kde=True)

    plt.axvline(1, color='red', linestyle='--', linewidth=2)

    plt.title(f'Distribution of Errors Out of Coverage Full (%)\n Dataset {ds}')
    plt.xlabel('Error %')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
    plt.savefig(
        rf"{prefix}\new_plots\errors_distribution_ds{ds}.png")
    plt.savefig(rf"{prefix}\new_plots\errors_distribution_ds{ds}.svg")

def plot_error_ratio_histogram(df, ds, bins=50):
    """
    Plots a histogram of ln(Full_Errors / Class_Errors)

    Parameters:
    - df: pandas DataFrame with 'Full_Errors' and 'Class_Errors' columns
    - bins: number of histogram bins (default: 50)
    """
    # Filter to avoid division by zero or negative values
    valid = (df['Full_Errors (%)'] > 0) & (df['Classified_Errors (%)'] > 0)
    ratios = np.log(df.loc[valid, 'Full_Errors (%)'] / df.loc[valid,
    'Classified_Errors (%)'])

    plt.figure(figsize=(8, 6))
    plt.hist(ratios, bins=bins, edgecolor='black')
    plt.xlabel('ln(Full_Errors (%) / Classified_Errors (%))')
    plt.ylabel('Frequency')
    plt.title(f"Histogram of ln(Full_Errors (%) / Classified_Errors (%)) \n "
              f"Dataset {ds}")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_full_only_coverage_histogram_dense(df, output_path="full_only_coverage_histogram.svg"):
    """
    Plots a normalized histogram of 'Full Only Coverage (%)' with 50 bins.
    Enlarged labels and ticks. Saves to SVG.
    """
    col = 'Full_Only_Coverage (%)'
    data = df[col].dropna()

    # Create 50 bins and normalize to percentages
    bins = np.linspace(0, data.max(), 70)
    weights = np.ones(len(data)) * 100.0 / len(data)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, edgecolor='black',
             color='gray', weights=weights)

    plt.xlabel('Full Only Coverage (%)', fontsize=16)
    plt.ylabel('Frequency (%)', fontsize=16)
    plt.title('Normalized Histogram (70 Bins)\nFull Only Coverage (%)',
              fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.savefig(output_path, format='svg')
    plt.show()


def plot_full_only_coverage_histogram_ranges(df):
    """
    Plots a normalized histogram of 'Full Only Coverage (%)' using 5 defined bins:
    0, (0-1], (1-5], (5-10], (10+]
    """
    col = 'Full_Only_Coverage (%)'
    data = df[col].dropna()

    # Define custom bins and labels
    bins = [-0.01, 0, 1, 5, 10, float('inf')]
    labels = ['0', '0-1', '1-5', '5-10', '10+']

    # Bin and normalize manually
    binned = pd.cut(data, bins=bins, labels=labels, right=True)
    counts = binned.value_counts(sort=False)
    percentages = counts / counts.sum() * 100

    # Plot
    plt.figure(figsize=(8, 6))
    percentages.plot(kind='bar', edgecolor='black', color='gray')

    plt.xlabel('Full Only Coverage (%)')
    plt.ylabel('Frequency (%)')
    plt.title(f'Normalized Histogram of Full Only Coverage (%)')
              # f'\nDataset {ds}')
    plt.grid(axis='y')
    plt.tight_layout()
    plt.show()



def plot_coverage_ranking(df, coverage_col='Full_Only_Coverage (%)'):
    """
    Plots rank vs. Full_Only_Coverage (%) values

    Parameters:
    - df: DataFrame with coverage column
    - coverage_col: name of the coverage column
    """
    # Drop missing values and sort
    values = df[coverage_col].dropna().sort_values(ascending=True).reset_index(drop=True)
    ranks = range(1, len(values) + 1)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(ranks, values, marker='o', linestyle='-', markersize=3)

    # Red dashed line at 5% coverage
    plt.axhline(5, color='red', linestyle='--', label='5% Coverage Threshold')

    # Set labels and title with larger font size
    plt.xlabel("Rank", fontsize=16)
    plt.ylabel("Full_Only_Coverage (%)", fontsize=16)
    plt.title("Ranking of Full_Only_Coverage (%) Across Organisms", fontsize=18)

    # Enlarge tick labels
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Legend with larger font size
    plt.legend(fontsize=14)

    plt.grid(True)
    plt.tight_layout()
    plt.savefig("coverage_ranking.svg", format='svg')
    plt.show()

def merge_datasets(ds_list, version=""):
    dfs = []
    for ds_sz in ds_list:
        csv_path = fr"{prefix}\new{version}_genome_coverage_comparison_{ds_sz[0]}_{ds_sz[1]}.csv"
        df = pd.read_csv(csv_path)
        df['ds'] = f"{ds_sz[0]}_{ds_sz[1]}"  # optional: add source identifier
        dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)
    return merged


def main():
    # # === Load the data ===
    # ds = 1
    # size = 128
    # csv_path = f"genome_coverage_comparison_{ds}_{size}.csv"
    # df = pd.read_csv(csv_path)
    caterogizes_to_plot = ["Coverage", "Mismatch", "Gap", "Errors"]
    ds_and_size = [
        # [1, 64],
        # [2, 64],
        # [3, 64],
                   [1, 128],
        [2, 128],
        [3, 128],
        [4, 128]
    ]
    for cat in caterogizes_to_plot:
        for ds_sz in ds_and_size:
            csv_path = fr"{prefix}\new_genome_coverage_comparison_{ds_sz[0]}_{ds_sz[1]}.csv"
            df = pd.read_csv(csv_path)
            df["Full_Mismatch_from_covered (%)"] = (df["Full_Mismatch (%)"] /
                                                    df["Full_Coverage (%)"]) * 100
            df["Classified_Mismatch_from_covered (%)"] = (df["Classified_Mismatch (%)"] /
                                                          df["Classified_Coverage (%)"]) * 100
            df["Full_Gap_from_covered (%)"] = (df["Full_Gap (%)"] / df["Full_Coverage (%)"]) * 100
            df["Classified_Gap_from_covered (%)"] = (df["Classified_Gap (%)"] /
                                                     df["Classified_Coverage (%)"]) * 100
            df["Full_Errors_from_covered (%)"] = (df["Full_Mismatch_from_covered (%)"]
                                                   + df["Full_Gap_from_covered (%)"])
            df["Classified_Errors_from_covered (%)"] = (
                        df["Classified_Mismatch_from_covered (%)"]
                        + df["Classified_Gap_from_covered (%)"])
            df.to_csv(csv_path)
            # scatter_plot2(cat, df, ds_sz[0], ds_sz[1])
            # box_plot(cat, df, ds_sz[0], ds_sz[1])
            # plot_errors_distrubution(df, ds_sz[0])

            # plot_error_ratio_histogram(df, ds_sz[0])
            # plot_full_only_coverage_histogram(df, ds_sz[0])

            # if cat == "Coverage":  # only call once per CSV
            #     plot_error_metrics(df, ds_sz[0], ds_sz[1], mode="of_Total")
            #     plot_error_metrics(df, ds_sz[0], ds_sz[1], mode="of_Genome")


    merged = merge_datasets(ds_and_size)
    plot_full_only_coverage_histogram_dense(merged)
    plot_full_only_coverage_histogram_ranges(merged)
    plot_coverage_ranking(merged)

    # for cat in caterogizes_to_plot:
        # csv_path = fr"{prefix}\new_genome_coverage_comparison_all.csv"
        # df = pd.read_csv(csv_path)
        # scatter_plot3(cat, merged)






if __name__ == "__main__":
    main()
