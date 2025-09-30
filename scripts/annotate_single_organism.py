import pandas as pd
import argparse
import os
from extract_reads_functionality_new_faster import (
    load_combined_organism_info,
    build_interval_tree_from_gtf,
    annotate_read
)

def annotate_for_organism(reads_csv, organism_name, bact_csv, euk_csv, output_dir, skip_existing=True):
    out_csv = os.path.join(output_dir, f"{organism_name}_annotated.csv")
    done_flag = os.path.join(output_dir, f"{organism_name}.done")

    if skip_existing and os.path.exists(done_flag):
        print(f"[SKIP] {organism_name} already completed.")
        return

    try:
        reads = pd.read_csv(reads_csv)
        reads = reads[reads["source"] == organism_name]
        if reads.empty:
            print(f"[SKIP] No reads found for {organism_name}")
            return

        gtf_map, tax_map = load_combined_organism_info(bact_csv, euk_csv)
        gtf_path = gtf_map.get(organism_name)
        tax_info = tax_map.get(organism_name, {})

        if not gtf_path or not os.path.exists(gtf_path):
            print(f"[FAIL] GTF not found for {organism_name}")
            return

        print(f"[START] Annotating {len(reads)} reads for {organism_name}")
        tree = build_interval_tree_from_gtf(gtf_path)
        annotated = []

        for _, row in reads.iterrows():
            ann = annotate_read(
                chrom=row["chromosome"],
                start=row["start"],
                end=row["end"],
                cds_info=row["cds_info"],
                interval_tree=tree
            )
            # Optionally rename field for clarity
            ann["cds_overlap_pct"] = ann.pop("overlap_percentage")

            annotated.append({**row.to_dict(), **tax_info, **ann})

        df_out = pd.DataFrame(annotated)
        df_out.to_csv(out_csv, index=False)

        with open(done_flag, "w") as f:
            f.write("DONE\n")

        print(f"[DONE] {organism_name}: saved {len(df_out)} reads to {out_csv}")

    except Exception as e:
        print(f"[ERROR] Failed for {organism_name}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, help="Path to the preprocessed reads CSV")
    parser.add_argument("--organism", required=True, help="Organism name to process")
    parser.add_argument("--bacteria_csv", required=True, help="CSV with bacteria metadata")
    parser.add_argument("--eukaryote_csv", required=True, help="CSV with eukaryote metadata")
    parser.add_argument("--output_dir", required=True, help="Directory to store per-organism annotated CSV")
    parser.add_argument("--skip_existing", action="store_true", help="Skip if output already exists")

    args = parser.parse_args()
    annotate_for_organism(
        reads_csv=args.reads,
        organism_name=args.organism,
        bact_csv=args.bacteria_csv,
        euk_csv=args.eukaryote_csv,
        output_dir=args.output_dir,
        skip_existing=args.skip_existing
    )
