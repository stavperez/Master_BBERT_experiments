import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from itertools import combinations

# === CONFIG ===
PARQUET_DIR = "data"
OUTPUT_DIR = "plot_pairs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

EUKS = {
    'Acanthochromis_polyacanthus': 'Acanthochromis_polyacanthus_R1_range_scores_len_emb.parquet',
    'Coffea_arabica': 'Coffea_arabica_R1_range_scores_len_emb.parquet',
    'Cyanidioschyzon_merolae': 'Cyanidioschyzon_merolae_strain_10D_R1_range_scores_len_emb.parquet',
    'Saccharomyces_paradoxus': 'Saccharomyces_paradoxus_R1_range_scores_len_emb.parquet',
    'Watersipora_subatra': 'Watersipora_subatra_R1_range_scores_len_emb.parquet',
    # 'Nakaseomyces_glabratus':
}

BACTS = {
    'Actinomadura_madurae': 'Actinomadura_madurae_R1_range_scores_len_emb.parquet',
    'Brevibacillus_brevis': 'Brevibacillus_brevis_R1_range_scores_len_emb.parquet',
    'Escherichia_coli': 'Escherichia_coli_R1_range_scores_len_emb.parquet',
    'Helicobacter_pylori': 'Helicobacter_pylori_R1_range_scores_len_emb.parquet',
    'Pseudomonas_aeruginosa': 'Pseudomonas_aeruginosa_R1_range_scores_len_emb.parquet',
}

ALL_ORGS = {**EUKS, **BACTS}

# Fixed pair color map
PAIR_COLORS = ['orange', 'green']  # Used when both are same domain
DOMAIN_COLORS = {'eukaryote': 'blue', 'bacteria': 'red'}  # Used when mixed

def get_domain(org):
    return 'eukaryote' if org in EUKS else 'bacteria'

# === FUNCTIONS ===
def flatten_embedding(embedding):
    return np.stack(embedding).astype(np.float32).flatten()

def load_embeddings(name, filename):
    path = os.path.join(PARQUET_DIR, filename)
    df = pd.read_parquet(path)
    df['flat_embedding'] = df['embedding'].apply(flatten_embedding)
    return df[['flat_embedding']].assign(label=name, domain=get_domain(name))

def plot_2d(df, x_col, y_col, title, filename_base):
    domains = df['domain'].unique()
    plt.figure(figsize=(10, 8))

    if len(domains) == 1:
        # same domain: assign different colors per label
        unique_labels = df['label'].unique()
        for i, label in enumerate(unique_labels):
            subset = df[df['label'] == label]
            plt.scatter(
                subset[x_col], subset[y_col],
                label=f"{label} ({subset['domain'].iloc[0]})",
                alpha=0.5,
                c=PAIR_COLORS[i % len(PAIR_COLORS)]
            )
    else:
        # mixed domains: use domain-based coloring
        for label in df['label'].unique():
            subset = df[df['label'] == label]
            domain = subset['domain'].iloc[0]
            plt.scatter(
                subset[x_col], subset[y_col],
                label=f"{label} ({domain})",
                alpha=0.5,
                c=DOMAIN_COLORS[domain]
            )

    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"pdf_{filename_base}.pdf"), dpi=600)
    # plt.savefig(os.path.join(OUTPUT_DIR, f"{filename_base}.png"), dpi=300)
    plt.close()

def main():
    org_pairs = list(combinations(ALL_ORGS.keys(), 2))
    for org1, org2 in org_pairs:
        try:
            df1 = load_embeddings(org1, ALL_ORGS[org1])
            df2 = load_embeddings(org2, ALL_ORGS[org2])
            df = pd.concat([df1, df2], ignore_index=True)

            embeddings = np.vstack(df['flat_embedding'].values)

            # PCA
            pca_result = PCA(n_components=2).fit_transform(embeddings)
            df['PCA1'], df['PCA2'] = pca_result[:, 0], pca_result[:, 1]
            plot_2d(df, 'PCA1', 'PCA2', f"PCA: {org1} vs {org2}", f"pca_{org1}_{org2}")

            # t-SNE
            tsne_result = TSNE(n_components=2, perplexity=30, n_iter=1000, random_state=42).fit_transform(embeddings)
            df['TSNE1'], df['TSNE2'] = tsne_result[:, 0], tsne_result[:, 1]
            plot_2d(df, 'TSNE1', 'TSNE2', f"t-SNE: {org1} vs {org2}", f"tsne_{org1}_{org2}")

            print(f"Finished: {org1} vs {org2}")

        except Exception as e:
            print(f"Failed: {org1} vs {org2} — {e}")


if __name__ == "__main__":
    main()
