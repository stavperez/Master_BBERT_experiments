import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

prefix = "unique_genus"
# prefix = "new"

def build_taxonomy_tree(df):
    """
    Builds a basic phylogenetic tree using taxonomic hierarchy:
    domain > phylum > class > order > family > genus > species

    Parameters:
    - df: DataFrame with taxonomic columns

    Returns:
    - ETE Tree object
    """
    # Define taxonomy ranks
    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus',
             'species']

    # Create a root tree node
    root = Tree(name="Root")

    # A mapping from taxonomic path to TreeNode
    node_dict = {}

    for _, row in df.iterrows():
        path = []
        current_node = root
        for rank in ranks:
            name = row.get(rank)
            if pd.isna(name):
                continue  # skip missing levels
            path.append(name)
            path_key = "|".join(path)
            if path_key not in node_dict:
                new_node = current_node.add_child(name=name)
                node_dict[path_key] = new_node
                current_node = new_node
            else:
                current_node = node_dict[path_key]

    return root



def build_taxonomy_tree_with_labels(df):
    """
    Builds a taxonomy tree and annotates internal nodes with their taxonomic rank
    (e.g., class, phylum, order).
    """
    ranks = ['domain', 'phylum', 'class', 'order',
             'family',
             # 'genus',
             'species'
             ]
    root = Tree(name="Root")
    node_dict = {}

    for _, row in df.iterrows():
        path = []
        current_node = root
        for rank in ranks:
            name = row.get(rank)
            if pd.isna(name):
                continue
            path.append(name)
            path_key = "|".join(path)
            if path_key not in node_dict:
                new_node = current_node.add_child(name=name)
                node_dict[path_key] = new_node
                current_node = new_node

                # Add label for key taxonomic ranks
                if rank in ['domain', 'phylum', 'class']: # , 'order', 'family'
                    tf = TextFace(f"{name} ({rank})", fsize=10,
                                  fgcolor='blue')
                    new_node.add_face(tf, column=0, position="branch-right")
            else:
                current_node = node_dict[path_key]

    return root


def show_tree(tree):
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = 'c'  # Circular tree
    tree.show(tree_style=ts)

def save_tree_to_pdf(tree, output_path="taxonomy_tree.pdf"):
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = 'c'  # or 'r' for rectangular
    tree.render(output_path, w=800, units='px', tree_style=ts)

def save_tree_to_newick(tree, output_path="tree1.nwk"):
    with open(output_path, "w") as f:
        f.write(tree.write(format=8))  # format=5 gives standard Newick


def highlight_high_coverage(tree, df, highlighted_col="Full_Only_Coverage (%)",
                            threshold=10):
    """
    Styles tree leaves with Full_Only_Coverage (%) > threshold

    Parameters:
    - tree: ETE3 tree object
    - df: DataFrame with species and Full_Only_Coverage (%)
    - threshold: highlight condition (default: 10)
    """
    # Create a mapping: species name → Full_Only_Coverage
    coverage_map = df.set_index('species')[highlighted_col].to_dict()

    for leaf in tree.iter_leaves():
        species_name = leaf.name
        coverage = coverage_map.get(species_name)
        if coverage is not None and coverage > threshold:
            ns = NodeStyle()
            ns["fgcolor"] = "red"
            ns["size"] = 20
            ns["shape"] = "sphere"
            # ns["bold"] = True
            leaf.set_style(ns)


from ete3 import NodeStyle, TextFace
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def show_colorbar(cmap, norm):
    """Displays a colorbar for the colormap used in tree highlighting."""
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)

    cb1 = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                       orientation='horizontal',
                       cax=ax)
    cb1.set_label('Full Only Coverage (%)')
    plt.show()

def highlight_full_only_coverage(tree, df,
                                 highlighted_col="Full_Only_Coverage (%)",
                                 cmap_name='coolwarm'):
    """
    Colors leaves by Full_Only_Coverage (%) using a color gradient,
    and writes the value next to each organism.

    Parameters:
    - tree: ETE3 tree object
    - df: DataFrame with 'species' and 'Full_Only_Coverage (%)'
    - cmap_name: matplotlib colormap (default: 'coolwarm')
    """
    # Build coverage map: species → value
    coverage_map = df.set_index('species')[highlighted_col].to_dict()
    values = list(coverage_map.values())

    # # Normalize to 0-100%
    # norm = mcolors.Normalize(vmin=0, vmax=50)
    # cmap = cm.get_cmap(cmap_name)

    # Custom colormap: green → yellow → red, clipped at 50
    cmap = mcolors.LinearSegmentedColormap.from_list("yellow_red", ["yellow", "red"])
    norm = mcolors.Normalize(vmin=0, vmax=50)  # clamp range to 0-50


    for leaf in tree.iter_leaves():
        species_name = leaf.name
        coverage = coverage_map.get(species_name)

        if coverage is not None:
            color = mcolors.to_hex(cmap(norm(coverage)))

            ns = NodeStyle()
            ns["fgcolor"] = color
            ns["size"] = 25
            ns["shape"] = "circle"
            leaf.set_style(ns)

            # Show actual numeric value
            face = TextFace(f"{coverage:.1f}%", fsize=9)
            leaf.add_face(face, column=1, position="aligned")

    # Optional: show colorbar as external plot
    show_colorbar(cmap, norm)



def save_itol_gradient(df, highlighted_col="Full_Only_Coverage (%)",
filepath="coverage_gradient.txt"):
    with open(filepath, "w") as f:
        f.write("DATASET_GRADIENT\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tFull Only Coverage\n")
        f.write("COLOR\t#ff0000\n")
        f.write(f"LEGEND_TITLE\t{highlighted_col}\n")
        f.write("LEGEND_SHAPES\t1\n")
        f.write("LEGEND_COLORS\t#ff0000\n")
        f.write("LEGEND_LABELS\t% coverage\n")
        f.write("FIELD_COLORS\t#ffff00,#ff0000\n")
        f.write("FIELD_LABELS\tLow,High\n")
        f.write("DATA\n")

        for _, row in df.iterrows():
            f.write(f"{row['species']}\t{row[highlighted_col]:.2f}\n")


def export_itol_branch_labels(tree, filepath="itol_branch_labels.txt"):
    """
    Exports branch labels for internal nodes (e.g., phylum) to iTOL format.
    Assumes phylum names appear as node names.
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("DATASET_BRANCHLABELS\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tPhylum labels\n")
        f.write("COLOR\t#000000\n")
        f.write("DATA\n")
        for node in tree.traverse():
            # Show only internal nodes with names at phylum level
            if not node.is_leaf() and hasattr(node, 'name') and node.name:
                # Optionally: check if it's a phylum based on depth or name structure
                if '(' in node.name and 'phylum' in node.name:
                    phylum = node.name.split(' ')[0]  # Get just the name part
                    f.write(f"{node.name}\t{phylum}\tbranch\n")

def export_phylum_branch_labels(tree, df, filepath="itol_phylum_labels.txt"):
    """
    Export a branch label (once per clade) for each phylum using iTOL's BRANCH_LABELS format.
    """
    # Map species to phylum
    species_to_phylum = df.set_index("species")["phylum"].dropna().to_dict()

    # Build phylum → species list mapping
    phylum_to_species = {}
    for species, phylum in species_to_phylum.items():
        phylum_to_species.setdefault(phylum, []).append(species)

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("DATASET_BRANCHLABELS\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tPhylum labels\n")
        f.write("COLOR\t#000000\n")
        f.write("DATA\n")

        for phylum, species_list in phylum_to_species.items():
            # Find leaves in the tree for this phylum
            try:
                leaves = [leaf for leaf in tree.iter_leaves() if leaf.name in species_list]
                if len(leaves) < 2:
                    continue  # skip singleton phyla
                mrca = tree.get_common_ancestor(leaves)
                f.write(f"{mrca.name or phylum}\t{phylum}\tbranch\n")
            except:
                print(f"Warning: Failed to find MRCA for phylum {phylum}")


def label_phylum_internal_nodes(tree, df):
    """
    Assign phylum names to internal nodes in the tree based on MRCA of species in each phylum.
    These names will appear in iTOL if 'show internal node names' is enabled.
    """
    species_to_phylum = df.set_index("species")["phylum"].dropna().to_dict()

    phylum_to_species = {}
    for species, phylum in species_to_phylum.items():
        phylum_to_species.setdefault(phylum, []).append(species)

    for phylum, species_list in phylum_to_species.items():
        leaves = [leaf for leaf in tree.iter_leaves() if leaf.name in species_list]
        if len(leaves) < 2:
            continue  # only label real clades
        try:
            mrca = tree.get_common_ancestor(leaves)
            # Only rename if it's unnamed
            if not mrca.name or mrca.name.startswith("node"):
                mrca.name = phylum
        except:
            print(f"Warning: Could not find MRCA for phylum {phylum}")

def save_itol_text_labels(df, highlighted_col="Full_Only_Coverage (%)",
                          filepath="coverage_labels.txt"):
    with open(filepath, "w", encoding="utf-8", newline='\n') as f:
        f.write("DATASET_TEXT\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tFull Only Coverage (%)\n")
        f.write("COLOR\t#000000\n")
        f.write("SHOW_INTERNAL\t0\n")
        f.write("DATA\n")
        for _, row in df.iterrows():
            species = str(row['species']).strip()
            value = row[highlighted_col]
            if pd.notna(value):
                f.write(f"{species}\t{value:.1f}\n")


    with open("coverage_labels.txt", "rb") as f:
        for _ in range(10):
            print(f.readline())


def error_phylogenetic_plots(ds_list, merged_selected, highlighted_col="Full_Only_Coverage (%)"):
    """
    Reads and merges multiple CSVs with the naming pattern:
    prefix\new1_genome_coverage_comparison_{ds[0]}_{ds[1]}.csv

    Parameters:
    - prefix (str): base path to the CSV files
    - ds_list (list of tuple): list of (ds0, ds1) tuples

    Returns:
    - merged DataFrame
    """
    # Read genome metadata
    genomes_csv_path = fr"{prefix}\genomes_bacteria_unique_genus.csv"
    genomes = pd.read_csv(genomes_csv_path)

    # Merge on Annotation ↔ organism_name
    joined = merged_selected.merge(genomes, left_on='organism_name', right_on='organism_name', how='left')
    joined = joined[joined["domain"] == "bacteria"]
    joined["FNR (%)"] = joined["FNR"].round(2).astype(str)
    joined = joined.sort_values(by="organism_name")
    joined.to_csv(fr"{highlighted_col}_all.csv",
                  index=False)
    # tree = build_taxonomy_tree(joined)
    tree = build_taxonomy_tree_with_labels(joined)
    # label_phylum_internal_nodes(tree, joined)
    # export_itol_branch_labels(tree, "itol_branch_labels.txt")
    # export_phylum_branch_labels(tree, joined, "itol_phylum_labels.txt")

    # highlight_high_coverage(tree, joined, threshold=10)
    highlight_full_only_coverage(tree, joined, highlighted_col)
    # plot_full_only_coverage_histogram_dense(joined)
    # plot_full_only_coverage_histogram_ranges(joined)
    save_tree_to_pdf(tree, f"phylum_labeled_tree_{highlighted_col}.pdf")
    save_tree_to_newick(tree, f"phylum_labeled_tree_{highlighted_col}.nwk")
    from itertools import islice

    print([leaf.name for leaf in islice(tree.iter_leaves(), 5)])

    save_itol_gradient(joined, highlighted_col, f"{highlighted_col}.txt")
    save_itol_text_labels(joined, highlighted_col, f"{highlighted_col}.txt")
    # show_tree(tree)


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
    caterogizes_to_plot = ["Coverage", "Mismatch", "Gap"]
    ds_and_size = [
        # [1, 64],
        # [2, 64],
        # [3, 64],
                   [1, 128],
        [2, 128],
        [3, 128],
        [4, 128]
    ]
    ######## For assembly results ##########
    merged = merge_datasets(ds_and_size)
    # Select relevant columns
    selected_cols = [
        'Annotation', 'Genome_Length', 'Full_Coverage (%)', 'Classified_Coverage (%)',
        'Union_Coverage (%)', 'Full_Only_Coverage (%)', 'Classified_Only_Coverage (%)', 'ds'
    ]
    merged_selected = merged[selected_cols]
    ####### For FN scores ########

    merged_selected = pd.read_csv(
        r"/assembly_outputs/organism_confusion_long.csv")
    merged_selected["FNR"] = (merged_selected["FN"] / (merged_selected["FN"]
                                                + merged_selected["TP"]) *
                                                      100)
    error_phylogenetic_plots(ds_and_size, merged_selected, "FNR")


if __name__ == "__main__":
    main()
