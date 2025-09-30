import pandas as pd
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import json
import os
import pickle


def clean_title_to_annotation(stitle):
    return stitle.split(",")[0].split()[0]  # May not be used now


def merged_len(tree):
    tree_copy = tree.copy()
    tree_copy.merge_overlaps()
    return sum(iv.end - iv.begin for iv in tree_copy)


def load_annotation_mappings(csv_path, ds, size):
    df = pd.read_csv(csv_path)
    start_row = (ds - 1) * size
    end_row = ds * size
    df = df.iloc[start_row:end_row]  # Select rows based on ds and size

    print("Annotation mapping file loaded:", len(df))

    sseqid_to_annotation = {}
    annotation_to_length = {}

    for _, row in df.iterrows():
        annotation = row["organism_name"]
        genome_length = row["genome_length"]
        sseqids = str(row["sequence_headers"]).split(";")
        sseqids = [s.split(".")[0] for s in sseqids]  # Normalize

        annotation_to_length[annotation] = genome_length
        for sseqid in sseqids:
            sseqid_to_annotation[sseqid] = annotation

    print("Total sseqid_to_annotation entries:", len(sseqid_to_annotation))
    return sseqid_to_annotation, annotation_to_length

def subtract_and_merge(tree1, tree2):
    result = IntervalTree()
    for iv in tree1:
        overlaps = tree2.overlap(iv.begin, iv.end)
        if not overlaps:
            result.add(iv)
        else:
            current = [iv]
            for ov in overlaps:
                next_current = []
                for seg in current:
                    if seg.end <= ov.begin or seg.begin >= ov.end:
                        next_current.append(seg)
                    else:
                        if seg.begin < ov.begin:
                            next_current.append(Interval(seg.begin, ov.begin))
                        if seg.end > ov.end:
                            next_current.append(Interval(ov.end, seg.end))
                current = next_current
            for seg in current:
                result.add(seg)
    result.merge_overlaps()
    return result

def build_trees_and_errors(df):
    trees = defaultdict(IntervalTree)
    mismatches = defaultdict(int)
    gaps = defaultdict(int)

    for _, row in df.iterrows():
        sseqid = str(row["sseqid"]).split(".")[0]
        start, end = sorted([int(row["sstart"]), int(row["send"])])
        trees[sseqid].add(Interval(start, end))
        mismatches[sseqid] += int(row["mismatch"])
        gaps[sseqid] += int(row["gapopen"])

    print("Trees built for:", len(trees))
    return trees, mismatches, gaps


def load_or_compute_pickle_triplet(prefix, df):
    paths = {
        "trees": f"{prefix}_trees.pkl",
        "mismatch": f"{prefix}_mismatch.pkl",
        "gaps": f"{prefix}_gaps.pkl"
    }

    if all(os.path.exists(p) for p in paths.values()):
        print(f"Loading all from pickle for prefix: {prefix}")
        trees = load_pickle(paths["trees"])
        mismatch = load_pickle(paths["mismatch"])
        gaps = load_pickle(paths["gaps"])
    else:
        print(f"Computing and saving pickles for prefix: {prefix}")
        trees, mismatch, gaps = build_trees_and_errors(df)
        save_pickle(trees, paths["trees"])
        save_pickle(mismatch, paths["mismatch"])
        save_pickle(gaps, paths["gaps"])

    return trees, mismatch, gaps


def save_pickle(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)


def load_pickle(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def compute_error_stats(df_error, annotation_to_sseqids, annotation_to_length, suffix="Full"):
    summary = {}
    df_error["sseqid"] = df_error["sseqid"].astype(str).str.split(".").str[0]

    for annot, sseqids in annotation_to_sseqids.items():
        subset = df_error[df_error["sseqid"].isin(sseqids)]
        se = subset["Single_Errors"].sum()
        ap = subset["Adjacent_Pairs"].sum()
        lr = subset["Longer_Runs"].sum()
        total_error = se + ap + lr
        genome_len = annotation_to_length.get(annot, 1)

        summary[annot] = {
            f"{suffix}_Single_Errors": int(se),
            f"{suffix}_Adjacent_Pairs": int(ap),
            f"{suffix}_Longer_Runs": int(lr),
            f"{suffix}_Total_Errors": int(total_error),
            f"{suffix}_%Single_Errors_of_Total": round(100 * se / total_error, 2) if total_error else None,
            f"{suffix}_%Adjacent_Pairs_of_Total": round(100 * ap / total_error, 2) if total_error else None,
            f"{suffix}_%Longer_Runs_of_Total": round(100 * lr / total_error, 2) if total_error else None,
            f"{suffix}_%Single_Errors_of_Genome": round(100 * se / genome_len, 4),
            f"{suffix}_%Adjacent_Pairs_of_Genome": round(100 * ap / genome_len, 4),
            f"{suffix}_%Longer_Runs_of_Genome": round(100 * lr / genome_len, 4),
        }

    return summary


def extract_error_intervals(row):
    qseq = row["qseq"]
    sseq = row["sseq"]
    sstart = int(row["sstart"])
    send = int(row["send"])

    # Determine direction
    step = 1 if send >= sstart else -1
    pos = sstart

    mismatch_intervals = []
    gap_intervals = []

    for q, s in zip(qseq, sseq):
        if s == "-":
            # gap in subject: skip position
            continue
        if q == "-":
            # insertion in subject (gap in query)
            gap_intervals.append(f"{pos}-{pos+1}")
        elif q != s:
            # mismatch
            mismatch_intervals.append(f"{pos}-{pos+1}")
        pos += step  # move along subject sequence

    return pd.Series([";".join(mismatch_intervals), ";".join(gap_intervals)])

def parse_intervals(series, sseqids):
    tree = IntervalTree()
    for _, row in series.iterrows():
        if row["sseqid"] not in sseqids:
            continue
        for entry in str(row["Mismatch_Intervals"]).split(";"):
            if "-" not in entry:
                continue
            start, end = map(int, entry.split("-"))
            tree.addi(start, end)
    return tree

def parse_gap_intervals(series, sseqids):
    tree = IntervalTree()
    for _, row in series.iterrows():
        if row["sseqid"] not in sseqids:
            continue
        for entry in str(row["Gap_Intervals"]).split(";"):
            if "-" not in entry:
                continue
            start, end = map(int, entry.split("-"))
            tree.addi(start, end)
    return tree

def tree_len(tree):
    tree_copy = tree.copy()
    tree_copy.merge_overlaps()
    return sum(iv.end - iv.begin for iv in tree_copy)

def adjusted_tree_len(tree, mismatch_dict, gap_dict, sseqids):
    tree_copy = tree.copy()
    tree_copy.merge_overlaps()
    total_len = sum(iv.end - iv.begin for iv in tree_copy)
    mismatches = sum(mismatch_dict.get(sid, 0) for sid in sseqids)
    gaps = sum(gap_dict.get(sid, 0) for sid in sseqids)
    return max(0, total_len - mismatches - gaps)


def main(ds, size, prefix="new", bact_file="genome_refseq_bacteria.csv"):
    # === File paths ===
    # classified_file = f"classified_bact_{ds}_{size}bact_{size}euk_blast_final.csv"
    # full_file = f"mix_{ds}_{size}bact_{size}euk_blast_final.csv"
    classified_file = fr"{prefix}\new_classified_bact_{ds}_{size}bact" \
                      f"_{size}euk_blast_with_error_types.csv"
    full_file = rf"{prefix}\mix_{ds}_{size}bact" \
                f"_{size}euk_blast_with_error_types.csv"

    genome_info_csv = f"{prefix}\{bact_file}"

    ###### Add Errors intervals
    input_file = f"{prefix}\mix_{ds}_128bact_128euk_blast_raw.tsv"
    output_file = f"{prefix}\mix_{ds}_128bact_128euk_blast_with_intervals.csv"
    if not os.path.exists(output_file):
    # if os.path.exists(output_file):

        df = pd.read_csv(input_file, sep="\t", header=None)
        df.columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq",
            "sseq"
        ]
        print("Processing mismatch and gap intervals...")
        df[["Mismatch_Intervals", "Gap_Intervals"]] = df.apply(
            extract_error_intervals, axis=1)
        # Drop sequences to reduce size
        df.drop(columns=["qseq", "sseq"], inplace=True)
        df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")
    #############

    # === Load data ===
    df_full = pd.read_csv(full_file)
    df_classified = pd.read_csv(classified_file)
    print("Loaded rows from full:", len(df_full))
    print("Loaded rows from classified:", len(df_classified))

    sseqid_to_annotation, annotation_to_length = load_annotation_mappings(
        genome_info_csv, int(ds), int(size))

    full_trees, full_mismatch, full_gaps = load_or_compute_pickle_triplet(
        rf"{prefix}\full_{ds}_{size}", df_full)
    class_trees, class_mismatch, class_gaps = load_or_compute_pickle_triplet(
        rf"{prefix}\new_class_{ds}_{size}", df_classified)

    annotation_to_sseqids = defaultdict(list)
    for sseqid, annot in sseqid_to_annotation.items():
        annotation_to_sseqids[annot].append(sseqid)

    print("Total annotations:", len(annotation_to_sseqids))
    results = []

    ####### error full only analysis
    # Load mismatch/gap interval data
    intervals_df = pd.read_csv(
        f"{prefix}\mix_{ds}_{size}bact_{size}euk_blast_with_intervals.csv")
    intervals_df["sseqid"] = \
    intervals_df["sseqid"].astype(str).str.split(".").str[0]  # normalize

    #######

    for annot, sseqids in annotation_to_sseqids.items():
        genome_length = annotation_to_length.get(annot)
        if not genome_length:
            continue

        def merge_trees(tree_dict, ids):
            merged = IntervalTree()
            for sid in ids:
                merged |= tree_dict.get(sid, IntervalTree())
            merged.merge_overlaps()
            return merged

        full_tree = merge_trees(full_trees, sseqids)
        class_tree = merge_trees(class_trees, sseqids)

        print(f"\nAnnotation: {annot}")
        print(f"  sseqids: {sseqids}")
        print(f"  full_tree length: {tree_len(full_tree)}")
        print(f"  class_tree length: {tree_len(class_tree)}")

        union_tree = (full_tree | class_tree).copy()
        union_tree.merge_overlaps()
        full_only_tree = subtract_and_merge(full_tree, class_tree)
        class_only_tree = subtract_and_merge(class_tree, full_tree)
        full_only_coverage = adjusted_tree_len(full_only_tree, full_mismatch,
                                        full_gaps, sseqids)
        classified_only_coverage = adjusted_tree_len(class_only_tree, class_mismatch,
                                        class_gaps, sseqids)
        full_mismatch_len = sum(full_mismatch.get(sid, 0) for sid in sseqids)
        full_gaps_len = sum(full_gaps.get(sid, 0) for sid in sseqids)
        class_mismatch_len = sum(class_mismatch.get(sid, 0) for sid in
                                 sseqids)
        class_gaps_len = sum(class_gaps.get(sid, 0) for sid in sseqids)

        result = {
            "Annotation": annot,
            "Genome_Length": genome_length,
            "Full_Coverage (%)": round(
                100 * adjusted_tree_len(full_tree, full_mismatch, full_gaps,
                                        sseqids) / genome_length, 2),
            "Classified_Coverage (%)": round(
                100 * adjusted_tree_len(class_tree, class_mismatch,
                                        class_gaps, sseqids) / genome_length,
                2),
            "Union_Coverage (%)": round(
                100 * adjusted_tree_len(union_tree, full_mismatch, full_gaps,
                                        sseqids) / genome_length, 2),
            "Full_Only_Coverage (%)": round(
                100 * full_only_coverage / genome_length,
                2),
            "Classified_Only_Coverage (%)": round(
                100 * classified_only_coverage / genome_length,
                2),
            "Full_Mismatch (%)": round(100 * full_mismatch_len /
                                       genome_length, 4),
            "Classified_Mismatch (%)": round(100 * class_mismatch_len / genome_length, 4),
            "Full_Gap (%)": round(100 * full_gaps_len / genome_length, 4),
            "Classified_Gap (%)": round(100 * class_gaps_len /
                                        genome_length, 4),
            "Full_Errors (%)": round(100 * (full_mismatch_len + full_gaps_len) /
                                       genome_length, 4),
            "Classified_Errors (%)": round(100 * (class_mismatch_len +
                                           class_gaps_len) / genome_length, 4),
        }
        results.append(result)

        ####### error full only analysis
        genome_rows = intervals_df[intervals_df["sseqid"].isin(sseqids)]
        mismatch_tree = parse_intervals(genome_rows, sseqids)
        gap_tree = parse_gap_intervals(genome_rows, sseqids)

        def total_overlap_len(tree1, tree2):
            return sum(
                min(iv1.end, iv2.end) - max(iv1.begin, iv2.begin)
                for iv1 in tree1
                for iv2 in tree2.overlap(iv1.begin, iv1.end)
                if min(iv1.end, iv2.end) > max(iv1.begin, iv2.begin)
            )

        full_only_mismatches = total_overlap_len(mismatch_tree,
                                                 full_only_tree)
        full_only_gaps = total_overlap_len(gap_tree, full_only_tree)
        class_only_mismatches = total_overlap_len(mismatch_tree,
                                                  class_only_tree)
        class_only_gaps = total_overlap_len(gap_tree, class_only_tree)

        # Use actual covered length (from Full_Coverage (%))
        full_covered_len = genome_length * result["Full_Coverage (%)"] / 100
        class_covered_len = genome_length * result["Classified_Coverage (%)"] / 100

        # Avoid divide-by-zero
        result["Full_only_mismatches_outof_full_cov"] = round((
                                                                full_only_mismatches /
                                                full_covered_len) * 100,
            4) if full_covered_len else 0
        result["Full_only_gaps_outof_full_cov"] = round((full_only_gaps /
                                              full_covered_len) * 100,
                                         4) if full_covered_len else 0
        result["Classified_only_mismatches_outof_class_cov"] = round((
                class_only_mismatches / class_covered_len) * 100,
            4) if class_covered_len else 0
        result["Classified_only_gaps_outof_class_cov"] = round((
                                class_only_gaps / class_covered_len) * 100,
            4) if class_covered_len else 0
        ###############

        ##### full only errors out of full only coverage / full errors
        # out of full coverage
        errors_outof_coverage_full_only = 0 if full_only_coverage == 0 else \
            (full_only_mismatches + full_only_gaps) / full_only_coverage
        result["errors_outof_coverage_full_only (%)"] = errors_outof_coverage_full_only
        errors_outof_coverage_full = 0 if full_covered_len == 0 else \
            (full_mismatch_len + full_gaps_len) / full_covered_len
        result["errors_outof_coverage_full (%)"] = errors_outof_coverage_full
        errors_outof_coverage_full_dist = 0 if errors_outof_coverage_full == 0 else \
            errors_outof_coverage_full_only / errors_outof_coverage_full
        result["errors_outof_coverage_full_dist (%)"] = errors_outof_coverage_full_dist

        errors_outof_coverage_classified_only = 0 if classified_only_coverage == 0 else \
            (class_only_mismatches + class_only_gaps) / classified_only_coverage
        result["errors_outof_coverage_classified_only (%)"] = errors_outof_coverage_classified_only
        errors_outof_coverage_classified = 0 if class_covered_len == 0 else \
            (class_mismatch_len + class_gaps_len) / class_covered_len
        result["errors_outof_coverage_classified (%)"] = errors_outof_coverage_classified
        errors_outof_coverage_classified_dist = 0 if errors_outof_coverage_classified == 0 else \
            errors_outof_coverage_classified_only / errors_outof_coverage_classified
        result["errors_outof_coverage_classified_dist (%)"] = errors_outof_coverage_classified_dist

    # === Additional error-type analysis ===
    error_file_full = f"{prefix}\mix_{ds}_{size}bact_{size}euk_blast_with_error_types.csv"
    error_file_class = f"{prefix}\classified_bact_{ds}_{size}bact_{size}euk_blast_with_error_types.csv"

    summary_full = {}
    summary_class = {}

    # if os.path.exists(error_file_full):
    #     df_error_full = pd.read_csv(error_file_full)
    #     print(f"Loaded full error-type data from: {error_file_full}")
    #     summary_full = compute_error_stats(df_error_full,
    #                                        annotation_to_sseqids,
    #                                        annotation_to_length,
    #                                        suffix="Full")
    #
    # else:
    #     print(f"Error-type file not found: {error_file_full}")
    #
    # if os.path.exists(error_file_class):
    #     df_error_class = pd.read_csv(error_file_class)
    #     print(f"Loaded classified error-type data from: {error_file_class}")
    #     summary_class = compute_error_stats(df_error_class,
    #                                         annotation_to_sseqids,
    #                                         annotation_to_length,
    #                                         suffix="Classified")
    #
    # else:
    #     print(f"Error-type file not found: {error_file_class}")

    # Add error stats to results
    for r in results:
        annot = r["Annotation"]
        r.update(summary_full.get(annot, {}))
        r.update(summary_class.get(annot, {}))

    print("\nTotal genomes processed:", len(results))

    output_df = pd.DataFrame(results)
    output_df.to_csv(fr"{prefix}\new1_genome_coverage_comparison_{ds}"
                     f"_{size}.csv",
                     index=False)
    print(fr"Saved: {prefix}\new1_genome_coverage_comparison_{ds}_{size}.csv")


if __name__ == "__main__":
    prefix = "unique_genus"
    bact_file = "genomes_bacteria_unique_genus.csv"
    ds_and_size = [
        # [1, 64],
        # # [2, 64],
        # [3, 64],
        #            [1, 128],
        # [2, 128],
        # [3, 128],
        [4, 128]
    ]
    for ds, size in ds_and_size:
        main(str(ds), str(size), prefix=prefix, bact_file=bact_file)

