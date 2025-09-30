import re
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time


def evaluate_classification(bact_file, euk_file):
    """
    Evaluate classification by comparing the true source label with the predicted classification.

    Args:
        bact_file (str): Path to the FASTQ file where reads were classified as bacteria.
        euk_file (str): Path to the FASTQ file where reads were classified as eukaryotes.

    Returns:
        dict: Counts of TP, FP, FN, TN.
    """
    conf_matrix = defaultdict(int)

    def process_fastq(file, predicted_label):
        """Helper function to process a FASTQ file and update the confusion matrix."""
        with open(file, "r") as f:
            for line in f:
                if line.startswith("@"):
                    source_match = re.search(r"Source:(\w+)", line)
                    if source_match:
                        source_label = source_match.group(1)
                        true_label = "Bacteria" if "Bact" in source_label else "Eukaryote"

                        if true_label == "Bacteria" and predicted_label == "Bacteria":
                            conf_matrix["TP"] += 1
                        elif true_label == "Eukaryote" and predicted_label == "Bacteria":
                            conf_matrix["FP"] += 1
                        elif true_label == "Bacteria" and predicted_label == "Eukaryote":
                            conf_matrix["FN"] += 1
                        elif true_label == "Eukaryote" and predicted_label == "Eukaryote":
                            conf_matrix["TN"] += 1

    # Process reads classified as bacteria
    process_fastq(bact_file, "Bacteria")

    # Process reads classified as eukaryote
    process_fastq(euk_file, "Eukaryote")

    return conf_matrix


def create_confusion_matrix(conf_matrix):
    """
    Create and display a confusion matrix.

    Args:
        conf_matrix (dict): Confusion matrix counts.

    Returns:
        pd.DataFrame: Confusion matrix as a DataFrame.
    """
    df = pd.DataFrame([[conf_matrix["TP"], conf_matrix["FP"]],
                       [conf_matrix["FN"], conf_matrix["TN"]]],
                      index=["Actual Bacteria", "Actual Eukaryote"],
                      columns=["Predicted Bacteria", "Predicted Eukaryote"])

    print("Confusion Matrix:")
    print(df)

    # Plot heatmap with green-to-red color scale
    plt.figure(figsize=(6, 4))
    sns.heatmap(df, annot=True, fmt="d", cmap="RdYlGn_r", center=0, linewidths=0)
    plt.title("Confusion Matrix")

    # Save the heatmap with a filename that contains the current time
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    filename = f"confusion_matrix_{timestamp}.png"
    plt.savefig(filename)
    print(f"Confusion matrix saved as {filename}")

    plt.show()

    return df


def extract_false_negatives(bact_file, euk_file, output_fn_file):
    """
    Extracts False Negative (FN) reads and writes them to a new FASTQ file.

    Args:
        bact_file (str): Path to the FASTQ file where reads were classified as bacteria.
        euk_file (str): Path to the FASTQ file where reads were classified as eukaryotes.
        output_fn_file (str): Path to output FASTQ file with FN reads.
    """
    fn_reads = []  # List to store FN reads

    with open(euk_file, "r") as f:
        lines = f.readlines()

    for i in range(0, len(lines),
                   4):  # Process FASTQ format in blocks of 4 lines
        header, sequence, plus, quality = lines[i:i + 4]

        # Extract source label
        source_match = re.search(r"Source:(\w+)", header)
        if source_match:
            source_label = source_match.group(1)
            true_label = "Bacteria" if "Bact" in source_label else "Eukaryote"

            # False Negative: True label is Bacteria, but classified as Eukaryote
            if true_label == "Bacteria":
                fn_reads.extend([header, sequence, plus, quality])

    # Write FN reads to a new FASTQ file
    with open(output_fn_file, "w") as fn_file:
        fn_file.writelines(fn_reads)

    print(f"False Negative (FN) reads saved to: {output_fn_file}")


def extract_binary_labels(bact_file, euk_file, error_type, output_file):
    """
    Extracts False Negative (FN) reads and writes them to a new FASTQ file.

    Args:
        bact_file (str): Path to the FASTQ file where reads were classified as bacteria.
        euk_file (str): Path to the FASTQ file where reads were classified as eukaryotes.
        error_type (str): in {"FN", "FP", "TP", "TN"}
        output_fn_file (str): Path to output FASTQ file with FN reads.
    """
    error_type_reads = []  # List to store FN reads
    file = bact_file if error_type in ["TP", "FP"] else euk_file
    with open(file, "r") as f:
        lines = f.readlines()

    for i in range(0, len(lines), 4):  # Process FASTQ format in blocks of
        header, sequence, plus, quality = lines[i:i + 4]

        # Extract source label
        source_match = re.search(r"Source:(\w+)", header)
        if source_match:
            source_label = source_match.group(1)
            true_label = "Bacteria" if "Bact" in source_label else "Eukaryote"

            if error_type in ["TP", "FN"] and true_label == "Bacteria":
                error_type_reads.extend([header, sequence, plus, quality])
            elif error_type in ["FP", "TN"] and true_label == "Eukaryote":
                error_type_reads.extend([header, sequence, plus, quality])

    with open(output_file, "w") as fn_file:
        fn_file.writelines(error_type_reads)

    print(f"False Negative (FN) reads saved to: {output_file}")


def extract_fn_scores(fn_fastq_file, euk_csv_file, output_csv_file):
    """
    Extracts False Negative read IDs from a FASTQ file, matches them with classification scores,
    and saves the results in a new CSV file.

    Args:
        fn_fastq_file (str): Path to the FASTQ file containing False Negative reads.
        euk_csv_file (str): Path to the CSV file containing classification scores.
        output_csv_file (str): Path to save the output CSV file.
    """

    # Load the Eukaryote classification results CSV
    euk_df = pd.read_csv(euk_csv_file)

    # Extract FN read IDs
    fn_read_ids = []
    with open(fn_fastq_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                read_id = line.split()[0][
                          1:]  # Extract ID (remove '@' character)
                fn_read_ids.append(read_id)

    # Convert FN read IDs to a DataFrame
    fn_df = pd.DataFrame(fn_read_ids, columns=["id_R1"])

    # Merge FN read IDs with classification scores
    merged_df = fn_df.merge(euk_df, on="id_R1", how="left")

    # Save to CSV
    merged_df.to_csv(output_csv_file, index=False)

    print(f"False Negative read scores saved to: {output_csv_file}")


def extract_scores(fastq_file, euk_score_file, bact_score_file, error_type, output_csv_file):
    """
    Extracts False Negative read IDs from a FASTQ file, matches them with classification scores,
    and saves the results in a new CSV file.

    Args:
        fn_fastq_file (str): Path to the FASTQ file containing False Negative reads.
        euk_csv_file (str): Path to the CSV file containing classification scores.
        error_type (str): in {"FN", "FP", "TP", "TN"}
        output_csv_file (str): Path to save the output CSV file.
    """

    # Load the Eukaryote classification results CSV
    score_file = euk_score_file if error_type in ["TN", "FN"] else \
        bact_score_file
    score_df = pd.read_csv(score_file)

    # Extract FN read IDs
    read_ids = []
    with open(fastq_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                read_id = line.split()[0][1:]
                read_ids.append(read_id)

    reads_df = pd.DataFrame(read_ids, columns=["id_R1"])

    # Merge FN read IDs with classification scores
    merged_df = reads_df.merge(score_df, on="id_R1", how="left")

    # Add column with error type
    merged_df["Category"] = error_type

    # Save to CSV
    merged_df.to_csv(output_csv_file, index=False)
    print(f"{error_type} scores saved to: {output_csv_file}")


def plot_kde(csv_file):
    """
    Reads the resulted CSV file, creates a KDE (smoothed histogram) of average_score,
    and adds a vertical line at x = 1.3386.

    Args:
        csv_file (str): Path to the CSV file containing matched FN reads and scores.
    """
    # Load the data
    df = pd.read_csv(csv_file)

    # Drop missing values
    df = df.dropna(subset=["average_score"])

    # Set the color palette
    sns.set_palette("Set2")
    set2_color = sns.color_palette("Set2")[
        0]  # Pick the first color from the Set2 palette

    # Plot KDE
    plt.figure(figsize=(8, 5))
    sns.kdeplot(df["average_score"], fill=True, color=set2_color, alpha=0.6)
    plt.axvline(x=1.3386, color='red', linestyle='dashed', linewidth=2,
                label='Threshold: 1.3386')

    # Labels and title
    plt.xlabel("Average Score")
    plt.ylabel("Density")
    plt.title("Smoothed KDE Distribution of Average Scores for FN Reads")
    plt.legend()

    plt.show()


def plot_kde_combined(df):
    """
    Creates a KDE plot for TP, TN, FP, and FN on the same graph using Set2 palette.

    Args:
        df (pd.DataFrame): DataFrame containing all categories.
    """
    # Drop missing values
    df = df.dropna(subset=["average_score"])

    # Set the color palette
    sns.set_palette("Set2")
    colors = sns.color_palette("Set2")

    # Plot KDE for each category
    plt.figure(figsize=(10, 6))
    categories = ["TP", "TN", "FP", "FN"]
    for i, category in enumerate(categories):
        subset = df[df["Category"] == category]
        sns.kdeplot(subset["average_score"], fill=True, color=colors[i],
                    alpha=0.6, label=category)

    # Add a reference line
    plt.axvline(x=1.3386, color='red', linestyle='dashed', linewidth=2,
                label='Threshold: 1.3386')
    plt.xlim(0.85, 1.60)  # Adjusts the x-axis to range from 0.85 to 1.60

    # Labels and title
    plt.xlabel("Average Score")
    plt.ylabel("Density")
    plt.title(
        "Smoothed KDE Distribution of Average Scores for TP, TN, FP, and FN Reads")
    plt.legend()

    plt.show()

def main():
    # File paths for classified reads
    source = r"C:\Users\stavp\CSE\Master\Lab\Project" \
             r"\genomes_100K_assembly_output\mixed_32_organisms"
    bact_file = source + r"\bact_mixed_32_organisms_R1.fastq"  # Reads
    # classified as bacteria
    euk_file = source + r"\euk_mixed_32_organisms_R1.fastq"  # Reads
    # classified as eukaryote

    # # Evaluate classification
    # conf_matrix = evaluate_classification(bact_file, euk_file)
    #
    # # Create and visualize confusion matrix
    # create_confusion_matrix(conf_matrix)

    # output_fn_file = source + r"\false_negative_reads.fastq"
    # extract_false_negatives(bact_file, euk_file, output_fn_file)

    euk_score_file = source + r"\euk_mixed_32_organisms.csv"
    bact_score_file = source + r"\bact_mixed_32_organisms.csv"

    # output_csv_file = source + r"\fn_reads_with_scores.csv"

    # extract_fn_scores(output_fn_file, euk_score_file, output_csv_file)
    # plot_kde(output_csv_file)

    all_error_type_csvs = []  # Store CSV file paths
    for error_type in ["TP", "FP", "FN", "TN"]:
        output_file = source + fr"\{error_type}_reads.fastq"
        # extract_binary_labels(bact_file, euk_file, error_type, output_file)
        output_csv_file = output_file.replace(".fastq", ".csv")

        extract_scores(output_file, euk_score_file, bact_score_file, error_type, output_csv_file)
        all_error_type_csvs.append(output_csv_file)

    all_error_type_dfs = [pd.read_csv(csv_file) for csv_file in all_error_type_csvs]
    all_error_type_df = pd.concat(all_error_type_dfs, ignore_index=True)

    # Save the final combined CSV
    combined_output_csv = source + r"\all_reads_with_scores.csv"
    all_error_type_df.to_csv(combined_output_csv, index=False)
    print(f"All categories combined and saved to: {combined_output_csv}")

    # Generate KDE plot for all categories
    plot_kde_combined(all_error_type_df)


def plot_confusion_matrix_new_data(ds):
    # Load the CSV file
    csv_path = fr'C:\Users\stavp\CSE\Master\Lab\Project\datasets\unique_genus\ds{ds}\all_reads_with_scores.csv'
    df = pd.read_csv(csv_path)

    # Count the occurrences of each category
    counts = df['Category'].value_counts().to_dict()

    # Ensure all categories are present
    for label in ['TP', 'TN', 'FP', 'FN']:
        counts.setdefault(label, 0)

    # Build confusion matrix:
    #            Predicted Positive | Predicted Negative
    # Actual Positive      TP       |        FN
    # Actual Negative      FP       |        TN

    conf_matrix = [
        [counts['TP'], counts['FN']],  # Actual Positive
        [counts['FP'], counts['TN']]  # Actual Negative
    ]
    print("TP: " + str(counts['TP']) + ", FN: " + str(counts['FN']) + ", "
                                                                      "FP: " +
          str(counts[
        'FP']) + ", TN: " + str(counts['TN']))

    # Plot
    plt.figure(figsize=(5, 4))
    sns.heatmap(conf_matrix, annot=True, fmt=',', cmap='YlOrRd',
                xticklabels=['Predicted Bacteria', 'Predicted Eukaryote'],
                yticklabels=['Actual Bacteria', 'Actual Eukaryote'])

    plt.title(f"Confusion Matrix for dataset {ds}")
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.tight_layout()
    plt.show()


def plot_confusion_matrix_from_results(dataset, conf_matrix):
    """
    This function takes the results of error_analysis_new.py from the
    cluster and plots it with enlarged fonts and saves as SVG.
    """
    # Plot
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        conf_matrix,
        annot=True,
        fmt=',',
        cmap='YlOrRd',
        xticklabels=['Predicted Bacteria', 'Predicted Eukaryote'],
        yticklabels=['Actual Bacteria', 'Actual Eukaryote'],
        annot_kws={"size": 16},
        cbar_kws={'shrink': 0.8}
    )

    plt.title(f"Confusion Matrix for {dataset}", fontsize=18)
    plt.xlabel("Predicted", fontsize=16)
    plt.ylabel("Actual", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()

    # Save as SVG
    plt.savefig(f"confusion_matrix_{dataset}.svg", format='svg')
    plt.show()
    plt.close()


if __name__ == "__main__":
    # ds = 1
    # for ds in [1, 2, 3, 4]:
    #     plot_confusion_matrix_new_data(ds)
    # plot_confusion_matrix_new_data(1)
    # main()
    """
    results from the cluster:
    """
    results_dict = \
        {"Dataset 1": {"TN": 4240640, "TP": 3570709,
                       "FN": 874296, "FP": 204903},
         "Dataset 2": {"TN": 4051386, "TP": 3563821,
                       "FN": 905667, "FP": 218526},
         "Dataset 3": {"TN": 4106729, "TP": 3514215, ####
                       "FN": 965800, "FP": 196278},####
         "Dataset 4": {"TN": 4161000, "TP": 3569229,
                       "FN": 910780, "FP": 211787}
         }
    total_results = {"TN": 0, "TP": 0, "FN": 0, "FP": 0}
    for dataset in results_dict.values():
        for key in total_results:
            total_results[key] += dataset[key]

    results_dict["All Datasets"] = total_results
    for result in results_dict.keys():
        plot_confusion_matrix_from_results(result, [
            [results_dict[result]["TP"], results_dict[result]["FN"]], # Actual Positive
            [results_dict[result]["FP"], results_dict[result]["TN"]] # Actual Negative
        ])

