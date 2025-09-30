import re
import pandas as pd

def load_source_label_maps(bact_csv, euk_csv):
    bact_df = pd.read_csv(bact_csv)
    euk_df = pd.read_csv(euk_csv)

    # Normalize column names if needed
    bact_sources = set(bact_df["organism_name"])
    euk_sources = set(euk_df["organism_name"])

    return bact_sources, euk_sources


def extract_id_to_source_map(fastq_files):
    id_to_source = {}

    for fastq_file in fastq_files:
        with open(fastq_file, "r") as f:
            for line in f:
                if line.startswith("@") and "Source:" in line:
                    read_id = line.split()[0].lstrip("@").split("/")[0]
                    source = line.strip().split()[-1].replace("Source:", "")
                    id_to_source[read_id] = source
    for i, (read_id, source) in enumerate(id_to_source.items()):
        print(f"{read_id} → {source}")
        if i == 9:
            break

    return id_to_source


def evaluate_using_source_lookup(bact_file, euk_file, all_scores_file, bact_source_file, euk_source_file, output_csv):
    # Load source classification sets
    bact_sources, euk_sources = load_source_label_maps(bact_source_file, euk_source_file)

    # Extract read_id -> source map
    id_to_source = extract_id_to_source_map([bact_file, euk_file])

    # Load predictions
    df = pd.read_csv(all_scores_file)

    sources = []
    true_labels = []
    confusion_types = []

    for _, row in df.iterrows():
        read_id = row["normalized_id"]
        is_bact_pred = row["is_bact"]

        source = id_to_source.get(read_id, None)
        if source is None:
            sources.append("Unknown")
            true_labels.append("Unknown")
            confusion_types.append("Unknown")
            continue

        sources.append(source)

        if source in bact_sources:
            true_label = "Bacteria"
        elif source in euk_sources:
            true_label = "Eukaryote"
        else:
            true_label = "Unknown"

        true_labels.append(true_label)

        pred_label = "Bacteria" if is_bact_pred else "Eukaryote"

        if true_label == "Bacteria" and pred_label == "Bacteria":
            confusion = "TP"
        elif true_label == "Eukaryote" and pred_label == "Bacteria":
            confusion = "FP"
        elif true_label == "Bacteria" and pred_label == "Eukaryote":
            confusion = "FN"
        elif true_label == "Eukaryote" and pred_label == "Eukaryote":
            confusion = "TN"
        else:
            confusion = "Unknown"

        confusion_types.append(confusion)

    df["source"] = sources
    df["true_label"] = true_labels
    df["confusion_type"] = confusion_types

    df.to_csv(output_csv, index=False)
    print(f"Annotated results saved to: {output_csv}")
    print("\nConfusion matrix counts:")
    print(df["confusion_type"].value_counts())

def main():
    ds = "3"
    base = f"/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/bbert_results{ds}_size128/"
    bact_file = base + "new_classification_output_bacteria_R1.fastq"
    euk_file = base + "new_classification_output_eukaryotes_R1.fastq"
    all_scores_file = base + "new_classification_output_merged_predictions.csv"
    output_csv = base + "annotated_classification_results.csv"

    bact_source_file = "/cs/usr/stavperez/sp/genomes_new/genomes_bacteria_unique_genus.csv"
    euk_source_file = "/cs/usr/stavperez/sp/genomes_new/genomes_eukaryotes_unique_genus.csv"

    evaluate_using_source_lookup(bact_file, euk_file, all_scores_file, bact_source_file, euk_source_file, output_csv)

if __name__ == "__main__":
    main()
