import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# === CONFIG ===

prefix = "unique_genus"
# ds = 1

GENOME_CSV = f"{prefix}/genomes_bacteria_unique_genus.csv"
#### unique genus ####
# ORGANISM = "Leclercia_adecarboxylata" # ds 1
# ORGANISM = "Salmonella_enterica_subsp._enterica_serovar_Gallinarum" # ds 1 fam
# ORGANISM = "Caldanaerobacter_subterraneus_subsp._tengcongensis_MB4" # ds 2
# ORGANISM = "Caldimonas_thermodepolymerans" # ds 2 fam
# ORGANISM = "Pauljensenia_hongkongensis" # ds 3
# ORGANISM = "Photorhabdus_laumondii_subsp._laumondii" # ds 3 fam
# ORGANISM = "Chloroflexus_aurantiacus_J-10-fl"  # ds 4
# ORGANISM = "Thermobifida_fusca" # ds 4 fam
WINDOW_START = 10_000
WINDOW_SIZE = 10_000
WINDOW_END = WINDOW_START + WINDOW_SIZE
OUTDIR = f"{prefix}/window_plots"
os.makedirs(OUTDIR, exist_ok=True)

#### new ####
# ORGANISM = "Caldilinea_aerophila_DSM_14535_=_NBRC_104270"
# ORGANISM = "Thermoflexus_hugenholtzii" # ds 2
# ORGANISM = "Caldicellulosiruptor_naganoensis" # ds 1

def plot_single_org_alignment_windows(ORGANISM, MIX_BLAST, CLASSIFIED_BLAST):
    print(f"=== Loading genome info for: {ORGANISM} ===")
    genome_df = pd.read_csv(GENOME_CSV)
    row = genome_df[genome_df["organism_name"] == ORGANISM]
    if row.empty:
        raise ValueError(f"Organism '{ORGANISM}' not found.")
    sseqids = [s.split(".")[0] for s in str(row.iloc[0]["sequence_headers"]).split(";") if s]
    print(f"Found {len(sseqids)} sseqids for organism.")

    # === Load and filter BLAST ===
    def load_and_filter(path, label):
        print(f"\nLoading BLAST file: {path}")
        cols = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq"
        ]
        df = pd.read_csv(path, sep="\t", names=cols)
        df["sseqid_short"] = df["sseqid"].str.split(".").str[0]
        df = df[df["sseqid_short"].isin(sseqids)].copy()
        print(f"{label} alignments loaded: {len(df)} rows after filtering for sseqids.")
        return df

    mix_df = load_and_filter(MIX_BLAST, "Mix")
    class_df = load_and_filter(CLASSIFIED_BLAST, "Classified")

    # === Normalize coordinates
    offset_map, offset = {}, 0
    print("\nComputing contig offsets:")
    for sid in sseqids:
        max_len = 0
        for df in [mix_df, class_df]:
            sub = df[df["sseqid_short"] == sid]
            if not sub.empty:
                max_len = max(max_len, sub[["sstart", "send"]].max().max())
        offset_map[sid] = offset
        print(f"  {sid}: offset={offset}, max_len={max_len}")
        offset += max_len + 1000

    def add_coords(df, label):
        print(f"\nAssigning global coordinates for: {label}")
        starts, ends = [], []
        for _, row in df.iterrows():
            off = offset_map[row["sseqid_short"]]
            start = min(row["sstart"], row["send"]) + off
            end = max(row["sstart"], row["send"]) + off
            starts.append(start)
            ends.append(end)
        df["global_start"] = starts
        df["global_end"] = ends
        print(f"{label} → coordinate range: {df['global_start'].min()} - {df['global_end'].max()}")
        return df

    mix_df = add_coords(mix_df, "Mix")
    class_df = add_coords(class_df, "Classified")

    # === Build color track
    def build_track(df, start, end, label):
        print(f"\nBuilding color track for {label} in window {start:,}-{end:,}")
        arr = np.full(end - start, "white", dtype=object)
        count = 0
        for _, row in df[(df["global_start"] < end) & (df["global_end"] > start)].iterrows():
            base = row["global_start"]
            for i, (q, s) in enumerate(zip(row["qseq"], row["sseq"])):
                gpos = base + i
                if gpos < start or gpos >= end:
                    continue
                local = gpos - start
                if q == '-' or s == '-':
                    arr[local] = "black"
                elif q != s:
                    if arr[local] != "black":
                        arr[local] = "red"
                else:
                    if arr[local] not in ("black", "red"):
                        arr[local] = "gray"
            count += 1
        print(f"  Processed {count} alignments for {label}")
        return arr

    mix_colors = build_track(mix_df, WINDOW_START, WINDOW_END, "Mix")
    class_colors = build_track(class_df, WINDOW_START, WINDOW_END, "Classified")

    # === Plot
    print("\nGenerating plot...")
    fig, ax = plt.subplots(figsize=(10, 2))
    color_map = {"white": "#ffffff", "gray": "#cccccc", "red": "#ff0000", "black": "#000000"}
    for i, (m, c) in enumerate(zip(mix_colors, class_colors)):
        ax.add_patch(Rectangle((WINDOW_START + i, 0.6), 1, 0.35, color=color_map[m]))
        ax.add_patch(Rectangle((WINDOW_START + i, 0.1), 1, 0.35, color=color_map[c]))

    ax.set_xlim(WINDOW_START, WINDOW_END)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.275, 0.775])
    ax.set_yticklabels(["Classified", "Mix"], fontsize=16)
    ax.set_title(f"{ORGANISM}\nPositions {WINDOW_START:,} to "
                 f"{WINDOW_END:,}", fontsize=20)
    ax.set_xlabel("Genome coordinate", fontsize=16)
    plt.tight_layout()
    outfile = f"{OUTDIR}/pdf_window_{ORGANISM}_{WINDOW_START}" \
              f"_{WINDOW_END}.png"
    # plt.savefig(outfile, dpi=150)
    # svg_outfile = outfile.replace(".png", ".svg")
    # plt.savefig(svg_outfile)
    pdf_outfile = outfile.replace(".png", ".pdf")
    plt.savefig(pdf_outfile)
    plt.close()
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    # ORGANISM = "Leclercia_adecarboxylata" # ds 1
    # ORGANISM = "Salmonella_enterica_subsp._enterica_serovar_Gallinarum" # ds 1 fam
    # ORGANISM = "Caldanaerobacter_subterraneus_subsp._tengcongensis_MB4" # ds 2
    # ORGANISM = "Caldimonas_thermodepolymerans" # ds 2 fam
    # ORGANISM = "Pauljensenia_hongkongensis" # ds 3
    # ORGANISM = "Photorhabdus_laumondii_subsp._laumondii" # ds 3 fam
    # ORGANISM = "Chloroflexus_aurantiacus_J-10-fl"  # ds 4
    # ORGANISM = "Thermobifida_fusca" # ds 4 fam
    ORGANSISM_PER_DS = {
        # "Leclercia_adecarboxylata": 1,
        #                 "Salmonella_enterica_subsp._enterica_serovar_Gallinarum": 1,
                        # "Caldanaerobacter_subterraneus_subsp._tengcongensis_MB4": 2,
                        "Caldimonas_thermodepolymerans": 2, # 330-360
                        # "Pauljensenia_hongkongensis": 3, # 90-100
                        # "Photorhabdus_laumondii_subsp._laumondii": 3, #90-100
                        # "Chloroflexus_aurantiacus_J-10-fl": 4,
                        # "Thermobifida_fusca": 4
    }
    for ORGANISM, ds in ORGANSISM_PER_DS.items():
        MIX_BLAST = f"{prefix}/mix_{ds}_128bact_128euk_blast_raw.tsv"
        CLASSIFIED_BLAST = f"{prefix}/classified_bact_" \
                           f"{ds}_128bact_128euk_blast_raw.tsv"
        # for window in range(0, 100_000, 10_000):
        for window in range(330000, 360000, 10000):

            WINDOW_START = window
            WINDOW_END = WINDOW_START + WINDOW_SIZE
            print(f"\n=== Processing window {WINDOW_START:,} to {WINDOW_END:,} ===")
            plot_single_org_alignment_windows(ORGANISM, MIX_BLAST,
                                              CLASSIFIED_BLAST)
