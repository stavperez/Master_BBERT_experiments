import pandas as pd
import argparse
import glob
import os

def normalize_read_id(read_id):
    return read_id.split("/")[0]

def aggregate_and_merge(input_dir, classification_csv, output_csv):
    # Step 1: Find all organism-specific annotated CSVs
    annotated_files = glob.glob(os.path.join(input_dir, "*_annotated.csv"))
    if not annotated_files:
        print("[ERROR] No annotated CSVs found in input directory.")
        return

    print(f"[INFO] Found {len(annotated_files)} annotated CSV files.")

    # Step 2: Load and combine all annotated CSVs
    dfs = []
    for f in annotated_files:
        df = pd.read_csv(f)
        dfs.append(df)
    df_combined = pd.concat(dfs, ignore_index=True)

    # Step 3: Normalize read IDs
    df_combined["normalized_id"] = df_combined["read_id"].apply(normalize_read_id)

    # Step 4: Load classification results
    df_class = pd.read_csv(classification_csv)
    print(f"[INFO] Loaded classification file with {len(df_class)} entries.")
    df_class.to_csv("/cs/usr/stavperez/sp/InSilicoSeq/dataset_unique_genus/output_reads_mix_256_orgs/annotated_reads_for_test_mix_1.csv", index=False)

    # Step 5: Merge on normalized read ID
    df_merged = pd.merge(df_combined, df_class, how="left", on="normalized_id")

    # Step 6: Drop the temporary column
    df_merged.drop(columns=["normalized_id"], inplace=True)

    # Step 7: Sort by read_idx to preserve FASTQ order
    if "read_idx" in df_merged.columns:
        df_merged.sort_values("read_idx", inplace=True)

    # Step 8: Save to file
    df_merged.to_csv(output_csv, index=False)
    print(f"[DONE] Saved merged output with {len(df_merged)} reads to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Directory with *_annotated.csv files")
    parser.add_argument("--classification_csv", required=True, help="Classification results CSV to merge")
    parser.add_argument("--output_csv", required=True, help="Final output CSV file path")

    args = parser.parse_args()

    aggregate_and_merge(
        input_dir=args.input_dir,
        classification_csv=args.classification_csv,
        output_csv=args.output_csv
    )
