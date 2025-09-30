import os
import pandas as pd
import gzip
from Bio import SeqIO
import argparse

def load_predictions(csv_r1, csv_r2, output_prefix):
    """Load predictions and probabilities from CSVs and compute final bacterial classification."""
    df_r1 = pd.read_csv(csv_r1, usecols=["id", "pred_bact", "bact_prob"])
    df_r2 = pd.read_csv(csv_r2, usecols=["id", "pred_bact", "bact_prob"])

    df_r1["normalized_id"] = df_r1["id"].str.replace(r"/1$", "", regex=True)
    df_r2["normalized_id"] = df_r2["id"].str.replace(r"/2$", "", regex=True)

    merged_df = pd.merge(df_r1, df_r2, on="normalized_id", suffixes=("_R1", "_R2"))

    def resolve_label(row):
        if row["pred_bact_R1"] == row["pred_bact_R2"]:
            return row["pred_bact_R1"] == 1
        else:
            avg_prob = (row["bact_prob_R1"] + row["bact_prob_R2"]) / 2
            return avg_prob > 0.5

    merged_df["is_bact"] = merged_df.apply(resolve_label, axis=1)

    output_csv_path = output_prefix + "_merged_predictions.csv"
    merged_df[["normalized_id", "is_bact"]].to_csv(output_csv_path, index=False)
    print(f"Saved merged prediction CSV to: {output_csv_path}")

    return merged_df.set_index("normalized_id")["is_bact"].to_dict()

def smart_open(file_path, mode="rt"):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)

def filter_fastq(input_r1, input_r2, predictions, output_prefix):
    output_r1_bact = output_prefix + "_bacteria_R1.fastq"
    output_r2_bact = output_prefix + "_bacteria_R2.fastq"
    output_r1_euk = output_prefix + "_eukaryotes_R1.fastq"
    output_r2_euk = output_prefix + "_eukaryotes_R2.fastq"

    print("Writing output FASTQ files to:")
    print(f"- {output_r1_bact}")
    print(f"- {output_r2_bact}")
    print(f"- {output_r1_euk}")
    print(f"- {output_r2_euk}")

    with smart_open(input_r1, "rt") as r1_handle, \
         smart_open(input_r2, "rt") as r2_handle, \
         open(output_r1_bact, "w") as r1_bact_out, \
         open(output_r2_bact, "w") as r2_bact_out, \
         open(output_r1_euk, "w") as r1_euk_out, \
         open(output_r2_euk, "w") as r2_euk_out:

        r1_iter = SeqIO.parse(r1_handle, "fastq")
        r2_iter = SeqIO.parse(r2_handle, "fastq")

        for record_r1, record_r2 in zip(r1_iter, r2_iter):
            norm_id = record_r1.id.split("/")[0]
            is_bact = predictions.get(norm_id)
            if is_bact is True:
                SeqIO.write(record_r1, r1_bact_out, "fastq")
                SeqIO.write(record_r2, r2_bact_out, "fastq")
            elif is_bact is False:
                SeqIO.write(record_r1, r1_euk_out, "fastq")
                SeqIO.write(record_r2, r2_euk_out, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Split FASTQ files based on pred_bact and bact_prob classification.")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ (gzipped)")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ (gzipped)")
    parser.add_argument("--csv_r1", required=True, help="CSV file with pred_bact and bact_prob for R1")
    parser.add_argument("--csv_r2", required=True, help="CSV file with pred_bact and bact_prob for R2")
    parser.add_argument("--output", required=True, help="Output prefix")
    args = parser.parse_args()

    print("Loading predictions...")
    predictions = load_predictions(args.csv_r1, args.csv_r2, args.output)
    print("Processing FASTQ files...")
    filter_fastq(args.r1, args.r2, predictions, args.output)
    print("Splitting complete.")

if __name__ == "__main__":
    main()
