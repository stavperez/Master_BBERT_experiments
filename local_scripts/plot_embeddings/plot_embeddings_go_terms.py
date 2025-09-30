import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# === CONFIG ===
INPUT_TYPE = '50_50_TP_FN1'
# INPUT_TYPE = 'reading_frame_-2'
PARQUET_FILES = {
    # 'Saccharomyces_paradoxus': 'data/Saccharomyces_paradoxus_R1_range_scores_len_emb.parquet',
    'Pseudomonas_aeruginosa': 'data/reads_by_confusion_type_R1_scores_len_emb.parquet',
    # 'Pseudomonas_aeruginosa': 'data/reading_frame_-2_scores_len_emb.parquet'
}
GO_SLIM_CSV = r'C:\Users\stavp\CSE\Master\Lab\Project\plot_embeddings\data' \
              r'\Pseudomonas_Saccharomyces_GO_info_with_frame.csv'
OUTPUT_DIR = "plots"
OUTPUT_JOINED_CSV = "embeddings_with_slim_labels.csv"

# Reduction settings
PERPLEXITY = 20
N_ITER = 1000
RANDOM_STATE = 42

# Plotting
ALPHA = 0.75
POINT_SIZE = 20
LABEL_FIELD = 'organism_label'
COLOR_FIELD_MULTI =  'GO_slim_categories' # 'CDS'
COLOR_FIELD_PRIMARY = 'GO_slim_primary'
COLOR_FIELD_CONFUSION = 'confusion_type'
LEGEND_MAX_CATS = 20
TITLE_PREFIX = 'Embeddings by GO Slim (Primary)'
TITLE_PREFIX_CONFUSION = 'Embeddings by Confusion Type'

os.makedirs(OUTPUT_DIR, exist_ok=True)

# === FUNCTIONS ===

def load_go_slim_csv(path):
    df = pd.read_csv(path)
    # df = df[['read_id', 'GO_slim_categories', 'confusion_type']].copy()
    df['GO_slim_categories'] = df['GO_slim_categories'].fillna('').astype(str)
    df['confusion_type'] = df['confusion_type'].fillna('Unknown')
    return df

def flatten_nested_embedding(emb):
    arr = np.stack(emb)
    return arr.flatten()

def load_and_prepare_embeddings(parquet_files):
    dfs = []
    for organism, path in parquet_files.items():
        df = pd.read_parquet(path)
        df['flat_embedding'] = df['embedding'].apply(flatten_nested_embedding)
        df[LABEL_FIELD] = organism
        df['read_id'] = df['id']
        dfs.append(df[['read_id', 'flat_embedding', LABEL_FIELD]])
    return pd.concat(dfs, ignore_index=True)

def compute_category_counts(series):
    counts = {}
    for s in series.dropna():
        for cat in [c.strip() for c in s.split(';') if c.strip()]:
            counts[cat] = counts.get(cat, 0) + 1
    return pd.Series(counts).sort_values(ascending=False)

def choose_primary_category(multi_str, global_counts):
    if not isinstance(multi_str, str) or not multi_str.strip():
        return np.nan
    cats = [c.strip() for c in multi_str.split(';') if c.strip()]
    if not cats:
        return np.nan
    return sorted(cats, key=lambda c: (len(c), -global_counts.get(c, 0)), reverse=True)[0]

def limit_to_top_categories(series, top_n):
    freq = series.value_counts(dropna=True)
    top = set(freq.head(top_n).index.tolist())
    mapping = {}
    for cat in freq.index:
        mapping[cat] = cat if cat in top else 'Other'
    return mapping

def apply_reduction(embeddings, method='pca'):
    if method == 'pca':
        reducer = PCA(n_components=2, random_state=RANDOM_STATE)
    elif method == 'tsne':
        reducer = TSNE(n_components=2, perplexity=PERPLEXITY, n_iter=N_ITER,
                       random_state=RANDOM_STATE, init='pca', learning_rate='auto')
    else:
        raise ValueError(f'Unknown method: {method}')
    return reducer.fit_transform(embeddings)

def generate_color_map(categories):
    strong_colors = [
        "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
        "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
        "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
        "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080"
    ]
    color_map = {}
    for i, cat in enumerate(sorted(categories)):
        if cat in ['Other', None, np.nan, 'nan']:
            color_map[cat] = '#d3d3d3'  # dark gray
        else:
            color_map[cat] = strong_colors[i % len(strong_colors)]
    color_map[np.nan] = '#d3d3d3'
    color_map['nan'] = '#d3d3d3'
    return color_map


def generate_confusion_color_map():
    return {
        'TP': '#1f77b4',     # blue
        'FP': '#ff7f0e',     # orange
        'TN': '#2ca02c',     # green
        'FN': '#d62728',     # red
        'Unknown': '#d3d3d3' # light gray
    }

def plot_2d(df, x_col, y_col, color_col, out_name):
    df_plot = df.copy()
    df_plot[color_col] = df_plot[color_col].fillna('nan')
    cats = sorted(df_plot[color_col].unique())
    color_map = generate_color_map(cats)

    plt.figure(figsize=(12, 10))
    for cat in cats:
        subset = df_plot[df_plot[color_col] == cat]
        plt.scatter(subset[x_col], subset[y_col],
                    s=POINT_SIZE, alpha=ALPHA,
                    c=[color_map[cat]], label=str(cat))

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f'{TITLE_PREFIX} - {x_col}/{y_col}')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
    plt.grid(True)
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, out_name)
    plt.savefig(out_path, dpi=300)
    print(f"Saved: {out_path}")
    plt.close()

def plot_2d_confusion(df, x_col, y_col, color_col, out_name):
    df_plot = df.copy()
    df_plot[color_col] = df_plot[color_col].fillna('Unknown')
    cats = sorted(df_plot[color_col].unique())
    color_map = generate_confusion_color_map()

    plt.figure(figsize=(10, 8))
    for cat in cats:
        subset = df_plot[df_plot[color_col] == cat]
        plt.scatter(subset[x_col], subset[y_col],
                    s=POINT_SIZE, alpha=ALPHA,
                    c=[color_map.get(cat, '#d3d3d3')], label=str(cat))

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f'{TITLE_PREFIX_CONFUSION} - {x_col}/{y_col}')
    plt.legend(title='Confusion Type', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, out_name)
    plt.savefig(out_path, dpi=300)
    print(f"Saved: {out_path}")
    plt.close()

# === MAIN ===

def main():
    print("Loading embeddings...")
    emb_df = load_and_prepare_embeddings(PARQUET_FILES)

    print("Loading GO Slim annotations...")
    slim_df = load_go_slim_csv(GO_SLIM_CSV)

    print("Merging on read_id...")
    df = emb_df.merge(slim_df, on='read_id', how='left')

    print("Computing GO Slim category counts...")
    counts = compute_category_counts(df[COLOR_FIELD_MULTI])
    print("Top 20 categories:")
    print(counts.head(20))

    print("Assigning primary category per read...")
    df[COLOR_FIELD_PRIMARY] = df[COLOR_FIELD_MULTI].apply(lambda x: choose_primary_category(x, counts))

    if LEGEND_MAX_CATS is not None:
        print(f"Limiting to top {LEGEND_MAX_CATS} categories...")
        mapping = limit_to_top_categories(df[COLOR_FIELD_PRIMARY], LEGEND_MAX_CATS)
        df[COLOR_FIELD_PRIMARY] = df[COLOR_FIELD_PRIMARY].map(mapping).fillna('Other')

    print("Saving merged CSV with slim labels...")
    df[['read_id', LABEL_FIELD, COLOR_FIELD_MULTI, COLOR_FIELD_PRIMARY, COLOR_FIELD_CONFUSION]].to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_JOINED_CSV), index=False)

    # print("Applying PCA...")
    emb_array = np.vstack(df['flat_embedding'].values)
    pca_result = apply_reduction(emb_array, method='pca')
    df['PCA1'], df['PCA2'] = pca_result[:, 0], pca_result[:, 1]
    plot_2d(df, 'PCA1', 'PCA2', COLOR_FIELD_PRIMARY, f'pca_{INPUT_TYPE}_by_slim.svg')
    plot_2d_confusion(df, 'PCA1', 'PCA2', COLOR_FIELD_CONFUSION,
                      f'pca_{INPUT_TYPE}_by_confusion.svg')

    print("Applying t-SNE...")
    tsne_result = apply_reduction(emb_array, method='tsne')
    df['TSNE1'], df['TSNE2'] = tsne_result[:, 0], tsne_result[:, 1]
    plot_2d(df, 'TSNE1', 'TSNE2', COLOR_FIELD_PRIMARY, f'tsne_'
                                                       f'{INPUT_TYPE}_by_slim.svg')
    plot_2d_confusion(df, 'TSNE1', 'TSNE2', COLOR_FIELD_CONFUSION,
                      f'tsne_{INPUT_TYPE}_by_confusion.svg')

    # Custom t-SNE plot with reading frame and confusion shape
    plot_tsne_by_reading_frame_and_confusion(df,
        f'tsne_{INPUT_TYPE}_by_reading_frame_and_confusion.svg')


    print("✅ Done.")


def plot_tsne_by_reading_frame_and_confusion(df, out_name):
    print("Plotting t-SNE by reading frame + confusion...")

    df_plot = df.copy()

    # Assign color groups
    def determine_color_group(row):
        if pd.notna(row.get('reading_frame')):
            return f"RF_{int(row['reading_frame'])}"
        elif row.get('CDS') == 'coding':
            return 'No_RF_Coding'
        elif row.get('CDS') == 'non-coding':
            return 'No_RF_NonCoding'
        else:
            return 'No_RF_Other'

    df_plot['color_group'] = df_plot.apply(determine_color_group, axis=1)

    # Color map
    color_groups = sorted(df_plot['color_group'].unique())
    color_map = generate_color_map(color_groups)

    # Shape map
    shape_map = {'TP': 'o', 'FN': 'x'}
    default_marker = '.'

    plt.figure(figsize=(12, 10))
    for color_group in color_groups:
        for confusion_type in ['TP', 'FN']:
            subset = df_plot[
                (df_plot['color_group'] == color_group) &
                (df_plot['confusion_type'] == confusion_type)
            ]
            if not subset.empty:
                plt.scatter(subset['TSNE1'], subset['TSNE2'],
                            s=POINT_SIZE, alpha=ALPHA,
                            c=[color_map[color_group]],
                            marker=shape_map.get(confusion_type, default_marker),
                            label=f"{color_group} ({confusion_type})")

    plt.xlabel('TSNE1')
    plt.ylabel('TSNE2')
    plt.title('t-SNE: Reading Frame + Confusion Type')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
    plt.grid(True)
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, out_name)
    plt.savefig(out_path, dpi=300)
    print(f"Saved: {out_path}")
    plt.close()

if __name__ == "__main__":
    main()
