import pandas as pd
import networkx as nx
from pyvis.network import Network
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
import re
from Bio import Phylo
import matplotlib.patches as mpatches

def generate_network(df, meta_df, threshold, output_file):
    G = nx.Graph()
    
    for sample in df.index:
        meta = meta_df[meta_df['sample_id'] == sample]
        
        if not meta.empty:
            m = meta.iloc[0]
            title_html = (
                f"<b>ID:</b> {sample}<br>"
                f"<b>Patient:</b> {m.get('patient_id', 'NA')}<br>"
                f"<b>Location:</b> {m.get('latitude', 'NA')}, {m.get('longitude', 'NA')}<br>"
                f"<b>Conclusion:</b> {m.get('conclusion', 'NA')}"
            )
            group = m.get('patient_id', 'Unknown') 
        else:
            title_html = sample
            group = 'Unknown'

        G.add_node(sample, title=title_html, group=group, label=sample)
    
    samples = df.index.tolist()
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            s1 = samples[i]
            s2 = samples[j]
            dist = df.iloc[i, j]
            
            if dist <= threshold:
                weight = 1.0 / (dist + 1)
                G.add_edge(s1, s2, weight=weight, title=f"{dist} SNPs", label=str(dist))

    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black")
    net.from_nx(G)
    
    net.set_options("""
    var options = {
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -8000,
          "springLength": 250,
          "springConstant": 0.04
        }
      }
    }
    """)
    
    net.write_html(output_file)

def generate_plots(df, output_prefix):
    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    distances = df.where(mask).stack().values
    
    # Histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(distances, bins=20, kde=True, color="skyblue")
    plt.title("Distribution of SNP Distances")
    plt.savefig(f"{output_prefix}_histogram.png")
    plt.close()
    
    # Heatmap
    n_samples = df.shape[0]
    fig_dim = max(12, n_samples * 0.6)
    plt.figure(figsize=(fig_dim, fig_dim))
    font_size = 10 if n_samples < 20 else 8
    
    sns.heatmap(
        df, 
        cmap="viridis", 
        annot=True, 
        fmt="d", 
        square=True,
        cbar_kws={"shrink": 0.8},
        annot_kws={"size": font_size}
    )
    
    plt.title("SNP Distance Heatmap")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap.png")
    plt.close()

    # Violin Plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(y=distances, color="lightblue", inner="quartile")
    plt.title("Plot of SNP Distances")
    plt.ylabel("SNP Distance")
    plt.savefig(f"{output_prefix}_violin.png")
    plt.close()

def get_lineage(conclusion):
    if pd.isna(conclusion): return "Unknown"
    match = re.search(r"Lineage lineage([\d\.]+)", str(conclusion))
    if match:
        return match.group(1)
    return "Unknown"

def generate_phylo_tree(tree_file, meta_df, output_file):
    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"Error reading tree file: {e}")
        return

    lineage_map = {}
    for _, row in meta_df.iterrows():
        l = get_lineage(row.get('conclusion'))
        main_clade = l.split('.')[0] if l != "Unknown" else "Unknown"
        lineage_map[row['sample_id']] = main_clade

    unique_clades = sorted(list(set(lineage_map.values())))
    palette = sns.color_palette("bright", len(unique_clades))
    color_map_mpl = {clade: color for clade, color in zip(unique_clades, palette)}
    def to_bio_color(mpl_color):
        return tuple(int(x * 255) for x in mpl_color)
        
    color_map_bio = {clade: to_bio_color(color) for clade, color in zip(unique_clades, palette)}
    
    gray_bio = (128, 128, 128)

    def color_clade(clade):
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_bio.get(l, gray_bio)
            clade.color = c
            return c
        else:
            child_colors = [color_clade(c) for c in clade]
            first_color = child_colors[0]
            if all(c == first_color for c in child_colors):
                clade.color = first_color
                return first_color
            else:
                clade.color = gray_bio 
                return gray_bio

    color_clade(tree.root)

    # Plot
    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    Phylo.draw(
        tree, 
        axes=ax, 
        do_show=False, 
        show_confidence=False,
        label_func=lambda x: x.name if x.is_terminal() else "",
        branch_labels=None
    )
    
    terminals = tree.get_terminals()
    for i, clade in enumerate(terminals):
        y_pos = i + 1
        x_pos = tree.distance(tree.root, clade)
        l = lineage_map.get(clade.name, "Unknown")
        c_mpl = color_map_mpl.get(l, "gray")
        ax.scatter(x_pos, y_pos, color=c_mpl, s=80, zorder=10, edgecolors='white', linewidth=0.5)

    handles = [mpatches.Patch(color=color_map_mpl[c], label=f"Lineage {c}") for c in unique_clades]
    plt.legend(handles=handles, title="TB Lineage", loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    
    plt.title("Phylogenetic Tree", fontsize=14)
    plt.xlabel("Genetic Distance", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help="Path to distance_matrix.tsv")
    parser.add_argument('--metadata', required=True, help="Path to metadata.tsv")
    parser.add_argument('--tree', required=False, help="Path to phylo_tree.nwk")
    parser.add_argument('--threshold', type=int, default=12, help="SNP threshold")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', index_col=0)
    df.index.name = "Sample"
    
    meta_df = pd.read_csv(args.metadata, sep='\t')

    generate_network(df, meta_df, args.threshold, "transmission_network.html")
    generate_plots(df, "stats")
    
    if args.tree:
        generate_phylo_tree(args.tree, meta_df, "phylo_tree_colored.png")

if __name__ == "__main__":
    main()
