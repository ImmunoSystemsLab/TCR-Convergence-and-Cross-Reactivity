import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

# =========================
# SECTION 1: TRBV usage table (β chain percentages)
# =========================
data_beta = {
    "TRBV6-4": [30.77, 1.09, 7.14, 0, 7.14],
    "TRBV10-3": [15.38, 1.44, 3.57, 0, 3.57],
    "TRBV30": [15.38, 1.04, 0, 0, 0],
    "TRBV19": [7.69, 3.36, 16.07, 18.18, 18.75],
    "TRBV27": [0, 9.69, 14.28, 18.18, 6.25],
    "TRBV10-2": [0, 0, 16.07, 0, 0]
}
epitopes = ["AAGIGILTV", "ELAGIGILTV", "IMDQVPFSV", "YLEPGPVTA", "KTWGQYWQV"]
df_beta = pd.DataFrame(data_beta, index=epitopes)

print("\n=== TRBV Usage Table (β chain, %) ===")
print(df_beta)

plt.figure(figsize=(8,5))
sns.heatmap(df_beta, annot=True, cmap="magma", cbar_kws={'label': 'Usage %'})
plt.title("TRBV Usage Heatmap (β chain)")
plt.show()

# =========================
# SECTION 2: TRAV usage table (α chain percentages)
# =========================
data_alpha = {
    "TRAV12-2": [46.15, 96.9, 10.71, 0, 0],
    "TRAV35": [15.38, 0, 0, 0, 0],
    "TRAV14": [0, 0, 7.14, 0, 12.5],
    "TRAV22": [0, 0, 5.36, 0, 0],
    "TRAV13-1": [0, 0, 0, 27.27, 0],
    "TRAV21": [0, 0, 0, 18.18, 18.75],
    "TRAV5": [0, 0, 0, 0, 25.0]
}
df_alpha = pd.DataFrame(data_alpha, index=epitopes)

print("\n=== TRAV Usage Table (α chain, %) ===")
print(df_alpha)

plt.figure(figsize=(8,5))
sns.heatmap(df_alpha, annot=True, cmap="viridis", cbar_kws={'label': 'Usage %'})
plt.title("TRAV Usage Heatmap (α chain)")
plt.show()

# =========================
# SECTION 3: Epitope–TRBV Network (β chain)
# =========================
epitope_trbv = {
    "AAGIGILTV": ["TRBV6-4", "TRBV10-3", "TRBV30", "TRBV19"],
    "ELAGIGILTV": ["TRBV27", "TRBV20-1", "TRBV12-3", "TRBV28", "TRBV7-9"],
    "IMDQVPFSV": ["TRBV19", "TRBV10-2", "TRBV27", "TRBV6-4"],
    "YLEPGPVTA": ["TRBV6-5", "TRBV7-3", "TRBV27", "TRBV19", "TRBV14", "TRBV12-4"],
    "KTWGQYWQV": ["TRBV12-3", "TRBV4-1", "TRBV20-1", "TRBV7-9", "TRBV13", "TRBV27",
                  "TRBV5-5", "TRBV6-5", "TRBV7-6", "TRBV16"]
}
print("\n=== Epitope–TRBV Connections (β chain) ===")
for epitope, genes in epitope_trbv.items():
    print(f"{epitope}: {', '.join(genes)}")

G = nx.Graph()
for epitope, genes in epitope_trbv.items():
    G.add_node(epitope, bipartite=0)
    for gene in genes:
        G.add_node(gene, bipartite=1)
        G.add_edge(epitope, gene)

plt.figure(figsize=(10,7))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=1500, font_size=10)
plt.title("Epitope–TRBV Network (β chain)")
plt.show()

# =========================
# SECTION 4: Epitope–TRAV Network (α chain)
# =========================
epitope_trav = {
    "AAGIGILTV": ["TRAV12-2", "TRAV35"],
    "ELAGIGILTV": ["TRAV12-2"],
    "IMDQVPFSV": ["TRAV12-2", "TRAV14", "TRAV22"],
    "YLEPGPVTA": ["TRAV13-1", "TRAV21", "TRAV30"],
    "KTWGQYWQV": ["TRAV5", "TRAV21", "TRAV14"]
}
print("\n=== Epitope–TRAV Connections (α chain) ===")
for epitope, genes in epitope_trav.items():
    print(f"{epitope}: {', '.join(genes)}")

G = nx.Graph()
for epitope, genes in epitope_trav.items():
    G.add_node(epitope, bipartite=0)
    for gene in genes:
        G.add_node(gene, bipartite=1)
        G.add_edge(epitope, gene)

plt.figure(figsize=(10,7))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw(G, pos, with_labels=True, node_color="lightgreen", node_size=1500, font_size=10)
plt.title("Epitope–TRAV Network (α chain)")
plt.show()

# =========================
# SECTION 5: Motif Overlap Analysis (β chain)
# =========================
epitope_motifs_beta = {
    "ELAGIGILTV": {"QYF","AFF","QFF","PQHF"},
    "AAGIGILTV": {"QYF","AFF","QHF","PQHF"},
    "IMDQVPFSV": {"QYF","AFF","QFF"},
    "YLEPGPVTA": {"QYF","AFF"},
    "KTWGQYWQV": {"QYF","LFF","QFF","PQHF"}
}
overlap_matrix_beta = pd.DataFrame(0, index=epitope_motifs_beta.keys(), columns=epitope_motifs_beta.keys())
for e1 in epitope_motifs_beta:
    for e2 in epitope_motifs_beta:
        if e1 != e2:
            overlap_matrix_beta.loc[e1, e2] = len(epitope_motifs_beta[e1].intersection(epitope_motifs_beta[e2]))

print("\n=== Motif Overlap Matrix (β chain) ===")
print(overlap_matrix_beta)

plt.figure(figsize=(8,6))
sns.heatmap(overlap_matrix_beta, annot=True, cmap="Blues", cbar_kws={'label': 'Shared Motifs'})
plt.title("Motif Overlap Heatmap (β chain)")
plt.show()

# =========================
# SECTION 6: Motif Overlap Analysis (α chain)
# =========================
epitope_motifs_alpha = {
    "ELAGIGILTV": {"IIF","LTF","LIF","LVF"},
    "AAGIGILTV": {"LIF","FYF","STF","LTF"},
    "IMDQVPFSV": {"LIF","LTF","FVF","LNF","QFF","QYF"},
    "YLEPGPVTA": {"QYF","LTF","AFF","LVF"},
    "KTWGQYWQV": {"LIF","LTF","IIF","MRF","PTF"}
}
overlap_matrix_alpha = pd.DataFrame(0, index=epitope_motifs_alpha.keys(), columns=epitope_motifs_alpha.keys())
for e1 in epitope_motifs_alpha:
    for e2 in epitope_motifs_alpha:
        if e1 != e2:
            overlap_matrix_alpha.loc[e1, e2] = len(epitope_motifs_alpha[e1].intersection(epitope_motifs_alpha[e2]))

print("\n=== Motif Overlap Matrix (α chain) ===")
print(overlap_matrix_alpha)

plt.figure(figsize=(8,6))






