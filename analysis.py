import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import Levenshtein
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

# =========================
# SECTION 1: بيانات ملخصة من النتائج
# =========================
tcr_counts = {
    "ELAGIGILTV": 2289,
    "AAGIGILTV": 13,
    "IMDQVPFSV": 56,
    "YLEPGPVTA": 11,
    "KTWGQYWQV": 16
}
df_counts = pd.DataFrame.from_dict(tcr_counts, orient="index", columns=["count"])

# =========================
# SECTION 2: Diversity metrics (كل الريبرتوار)
# =========================
p = df_counts["count"] / df_counts["count"].sum()
shannon_entropy = -np.sum(p * np.log2(p))
simpson_index = 1 - np.sum(p**2)

print("\n=== Diversity Metrics (global) ===")
print(f"Shannon diversity index: {shannon_entropy:.3f}")
print(f"Simpson diversity index: {simpson_index:.3f}")

# =========================
# SECTION 3: TRBV enrichment (β chain)
# =========================
trbv_usage = {
    "ELAGIGILTV": {"TRBV27":222, "TRBV20-1":165, "TRBV12-3":118},
    "AAGIGILTV": {"TRBV6-4":4, "TRBV10-3":2, "TRBV30":2},
    "IMDQVPFSV": {"TRBV19":9, "TRBV10-2":9, "TRBV27":8},
    "YLEPGPVTA": {"TRBV6-5":3, "TRBV7-3":2, "TRBV27":2},
    "KTWGQYWQV": {"TRBV12-3":3, "TRBV4-1":3, "TRBV20-1":2}
}
genes = set().union(*[usage.keys() for usage in trbv_usage.values()])
table = []
for gene in genes:
    row = [trbv_usage[e].get(gene,0) for e in trbv_usage]
    table.append(row)

chi2, pval, dof, expected = chi2_contingency(table)
print(f"\nChi-square test for TRBV enrichment (β chain): p={pval:.4e}")

# =========================
# SECTION 4: Motif counts (β chain)
# =========================
motif_counts_beta = {
    "ELAGIGILTV": {"QYF": 792, "AFF": 322, "PQHF": 220, "LFF": 0},
    "AAGIGILTV": {"QYF": 3, "AFF": 4, "PQHF": 2, "LFF": 0},
    "IMDQVPFSV": {"QYF": 19, "AFF": 20, "PQHF": 0, "LFF": 0},
    "YLEPGPVTA": {"QYF": 9, "AFF": 2, "PQHF": 0, "LFF": 0},
    "KTWGQYWQV": {"QYF": 4, "AFF": 0, "PQHF": 2, "LFF": 4}
}
motifs = ["QYF","AFF","PQHF","LFF"]
epitopes = list(motif_counts_beta.keys())
table_beta = []
for motif in motifs:
    row = [motif_counts_beta[e].get(motif,0) for e in epitopes]
    table_beta.append(row)
df_beta_table = pd.DataFrame(table_beta, index=motifs, columns=epitopes)
print("\n=== Motif counts per epitope (β chain) ===")
print(df_beta_table)
chi2, pval, dof, expected = chi2_contingency(df_beta_table)
print(f"\nChi-square test for motif enrichment (β chain): p={pval:.4e}")














# =========================
# SECTION 2: TRAV usage (α chain)
# =========================
trav_usage = {
    "ELAGIGILTV": {"TRAV12-2":2218},
    "AAGIGILTV": {"TRAV12-2":6, "TRAV35":2},
    "IMDQVPFSV": {"TRAV12-2":6, "TRAV14":4, "TRAV22":3},
    "YLEPGPVTA": {"TRAV13-1":3, "TRAV21":2, "TRAV30":2},
    "KTWGQYWQV": {"TRAV5":4, "TRAV21":3, "TRAV14":2}
}

genes_alpha = list(set().union(*[usage.keys() for usage in trav_usage.values()]))

table_alpha = []
for gene in genes_alpha:
    row = [trav_usage[e].get(gene,0) for e in epitopes]
    table_alpha.append(row)

df_trav = pd.DataFrame(table_alpha, index=genes_alpha, columns=epitopes)
print("\n=== TRAV usage (α chain) ===")
print(df_trav)

chi2, pval, dof, expected = chi2_contingency(df_trav)
print(f"\nChi-square test for TRAV enrichment (α chain): p={pval:.4e}")


# =========================
# SECTION 3: TRAJ usage (α chain)
# =========================
traj_usage = {
    "ELAGIGILTV": {"TRAJ3":1875},
    "AAGIGILTV": {"TRAJ23":4, "TRAJ49":2},
    "IMDQVPFSV": {"TRAJ42":4, "TRAJ26":3},
    "YLEPGPVTA": {"TRAJ52":5, "TRAJ53":2},
    "KTWGQYWQV": {"TRAJ34":2, "TRAJ42":2, "TRAJ30":2}
}

genes_traj = list(set().union(*[usage.keys() for usage in traj_usage.values()]))

table_traj = []
for gene in genes_traj:
    row = [traj_usage[e].get(gene,0) for e in epitopes]
    table_traj.append(row)

df_traj = pd.DataFrame(table_traj, index=genes_traj, columns=epitopes)
print("\n=== TRAJ usage (α chain) ===")
print(df_traj)

chi2, pval, dof, expected = chi2_contingency(df_traj)
print(f"\nChi-square test for TRAJ enrichment (α chain): p={pval:.4e}")

# =========================
# SECTION 5: Motif counts (α chain)
# =========================
alpha_motif_counts = {
    "ELAGIGILTV": {"IIF": 1883, "LTF": 155, "LIF": 62, "LVF": 37},
    "AAGIGILTV": {"LIF": 4, "FYF": 2, "STF": 1, "LTF": 1},
    "IMDQVPFSV": {"LIF": 6, "LTF": 4, "FVF": 3, "LNF": 3, "QFF": 2, "QYF": 2},
    "YLEPGPVTA": {"QYF": 5, "LTF": 3, "AFF": 2, "LVF": 1},
    "KTWGQYWQV": {"LIF": 6, "LTF": 3, "IIF": 2, "MRF": 1, "PTF": 1}
}

# =========================
# SECTION 2: بناء جدول contingency
# =========================
motifs_alpha = list(set().union(*[usage.keys() for usage in alpha_motif_counts.values()]))  # حولناها لـlist
epitopes = list(alpha_motif_counts.keys())

table_alpha = []
for motif in motifs_alpha:
    row = [alpha_motif_counts[e].get(motif,0) for e in epitopes]
    table_alpha.append(row)

df_alpha_table = pd.DataFrame(table_alpha, index=motifs_alpha, columns=epitopes)
print("\n=== Motif counts per epitope (α chain) ===")
print(df_alpha_table)

# =========================
# SECTION 3: Chi-square test
# =========================
chi2, pval, dof, expected = chi2_contingency(df_alpha_table)
print(f"\nChi-square test for α motif enrichment across epitopes: p={pval:.4e}")

# =========================
# SECTION 6: CDR3 similarity clustering (Levenshtein motifs)
# =========================
cdr3_seqs = [
    "QYF","EQYF","TQYF","AFF","EAFF","QFF","EQFF","YEQYF","TEAFF","PQHF",
    "GTEAFF","FGTEAFF","DTQYF","LFF","ELFF","GELFF"
]
dist_matrix = np.zeros((len(cdr3_seqs), len(cdr3_seqs)))
for i, s1 in enumerate(cdr3_seqs):
    for j, s2 in enumerate(cdr3_seqs):
        dist_matrix[i,j] = Levenshtein.distance(s1, s2)
condensed_dist = squareform(dist_matrix)
linkage_matrix = linkage(condensed_dist, method="average")
plt.figure(figsize=(12,6))
dendrogram(linkage_matrix, labels=cdr3_seqs, leaf_rotation=90)
plt.title("CDR3 motif similarity clustering (Levenshtein distance)")
plt.show()

clusters = fcluster(linkage_matrix, t=3, criterion="distance")
cluster_df = pd.DataFrame({"CDR3": cdr3_seqs, "Cluster": clusters})
print("\n=== CDR3 Clusters (threshold=3) ===")
print(cluster_df.sort_values("Cluster"))
