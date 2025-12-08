import pandas as pd
import networkx as nx
from pyvis.network import Network
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np

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
    plt.figure(figsize=(12, 10))
    sns.heatmap(df, cmap="viridis", annot=True, fmt="d")
    plt.title("SNP Distance Heatmap")
    plt.savefig(f"{output_prefix}_heatmap.png")
    plt.close()

    # Violin Plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(y=distances, color="lightblue", inner="quartile")
    plt.title("Plot of SNP Distances")
    plt.ylabel("SNP Distance")
    plt.savefig(f"{output_prefix}_violin.png")
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help="Path to distance_matrix.tsv")
    parser.add_argument('--metadata', required=True, help="Path to metadata.tsv")
    parser.add_argument('--threshold', type=int, default=12, help="SNP threshold")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', index_col=0)
    df.index.name = "Sample"
    
    meta_df = pd.read_csv(args.metadata, sep='\t')

    generate_network(df, meta_df, args.threshold, "transmission_network.html")
    generate_plots(df, "stats")

if __name__ == "__main__":
    main()
