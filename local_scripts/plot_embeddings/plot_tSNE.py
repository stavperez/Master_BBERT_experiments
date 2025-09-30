import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
# import umap
from matplotlib.gridspec import GridSpec

# === CONFIG ===
PARQUET_DIR = "data"
PARQUET_FILES = {
    # 'Saccharomyces_paradoxus': 'data\Saccharomyces_paradoxus_R1_range_scores_len_emb'
    #                '.parquet',
    'Pseudomonas_aeruginosa': 'data\Pseudomonas_aeruginosa_R1_range_scores_len_emb'
                    '.parquet',
}

SINGLE_FILES = {
    'Acanthochromis_polyacanthus': 'Acanthochromis_polyacanthus_R1_range_scores_len_emb.parquet',
    'Coffea_arabica': 'Coffea_arabica_R1_range_scores_len_emb.parquet',
    'Cyanidioschyzon_merolae': 'Cyanidioschyzon_merolae_strain_10D_R1_range_scores_len_emb.parquet',
    'Saccharomyces_paradoxus': 'Saccharomyces_paradoxus_R1_range_scores_len_emb.parquet',
    'Watersipora_subatra': 'Watersipora_subatra_R1_range_scores_len_emb.parquet',
    # 'Nakaseomyces_glabratus':
    'Actinomadura_madurae': 'Actinomadura_madurae_R1_range_scores_len_emb.parquet',
    'Brevibacillus_brevis': 'Brevibacillus_brevis_R1_range_scores_len_emb.parquet',
    'Escherichia_coli': 'Escherichia_coli_R1_range_scores_len_emb.parquet',
    'Helicobacter_pylori': 'Helicobacter_pylori_R1_range_scores_len_emb.parquet',
    'Pseudomonas_aeruginosa': 'Pseudomonas_aeruginosa_R1_range_scores_len_emb.parquet',
}

OUTPUT_DIR = "plots"
LABEL_FIELD = "label"
PERPLEXITY = 20
N_ITER = 1000

# === COLOR CONFIG ===
colormap = sns.color_palette("Set1", 10)
COLORS = {
    'S.paradoxus': colormap[0],
    'Pseudomonas_aeruginosa': colormap[2],
}
ALPHA = 0.4

# === FUNCTIONS ===
def load_and_prepare_data(parquet_files):
    all_dfs = []
    for label, path in parquet_files.items():
        df = pd.read_parquet(path)
        df[LABEL_FIELD] = label

        def flatten_nested_embedding(emb):
            arr = np.stack(emb)  # shape: (102, 768)
            return arr.flatten()  # shape: (78336,)

        df['flat_embedding'] = df['embedding'].apply(flatten_nested_embedding)

        all_dfs.append(df)
    return pd.concat(all_dfs, ignore_index=True)

def load_single_data(label, path):
    path = os.path.join(PARQUET_DIR, path)
    df = pd.read_parquet(path)
    df[LABEL_FIELD] = label

    def flatten_nested_embedding(emb):
        arr = np.stack(emb)  # shape: (102, 768)
        return arr.flatten()  # shape: (78336,)

    df['flat_embedding'] = df['embedding'].apply(flatten_nested_embedding)

    return df
def apply_reduction(embeddings, method="pca"):
    if method == "pca":
        reducer = PCA(n_components=2)
    elif method == "tsne":
        reducer = TSNE(n_components=2, perplexity=PERPLEXITY, n_iter=N_ITER, random_state=42)
    elif method == "umap":
        reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    else:
        raise ValueError(f"Unknown method: {method}")
    return reducer.fit_transform(embeddings)

def plot_2d(df, x_col, y_col, title, output_file):
    plt.figure(figsize=(10, 8))
    plt.scatter(df[x_col].values, df[y_col].values, alpha=ALPHA,
                c="grey")  #colormap[0])
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, output_file), dpi=600, format='pdf')
    print(f"Saved: {output_file}")
    plt.close()

def main(label, path):
    df = load_single_data(label, path)
    # df = load_and_prepare_data(PARQUET_FILES)

    # embeddings = np.vstack(df['embedding'].values)
    embeddings = np.vstack(
        df['flat_embedding'].values)  # shape: (num_samples, 78336)
    print(df['flat_embedding'].iloc[0].shape)  # should print (78336,)
    print(embeddings.shape)  # should print (N, 78336)

    # PCA
    pca_result = apply_reduction(embeddings, method='pca')
    df['PCA1'], df['PCA2'] = pca_result[:, 0], pca_result[:, 1]
    plot_2d(df, 'PCA1', 'PCA2', f"PCA: {label}", f"pca_{label}.pdf")

    # t-SNE
    tsne_result = apply_reduction(embeddings, method='tsne')
    df['TSNE1'], df['TSNE2'] = tsne_result[:, 0], tsne_result[:, 1]
    plot_2d(df, 'TSNE1', 'TSNE2', f"t-SNE: {label}", f"tsne_{label}.pdf")

    # UMAP
    # umap_result = apply_reduction(embeddings, method='umap')
    # df['UMAP1'], df['UMAP2'] = umap_result[:, 0], umap_result[:, 1]
    # plot_2d(df, 'UMAP1', 'UMAP2', "UMAP of Embeddings", "umap_plot.svg")


if __name__ == "__main__":
    for label, path in SINGLE_FILES.items():
        main(label, path)
