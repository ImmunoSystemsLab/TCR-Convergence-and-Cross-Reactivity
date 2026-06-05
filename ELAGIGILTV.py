import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

# Load dataset
df = pd.read_csv("VDJdb_ELAGIGILTV_TRA_TRB_Input.csv", sep="\t")
print("Total TCRs:", df.shape[0])

# =========================
# SECTION A: β-chain analysis
# =========================
if "CDR3b" in df.columns:
    df["cdr3b_length"] = df["CDR3b"].str.len()
    print("\nDistribution of CDR3β lengths:", df["cdr3b_length"])

    plt.figure(figsize=(8,5))
    sns.histplot(df["cdr3b_length"], bins=15, kde=True, color="red")
    plt.title("CDR3β length distribution")
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.show()

if "TRBV" in df.columns:
    v_usage_b = df["TRBV"].value_counts()
    plt.figure(figsize=(10,5))
    sns.barplot(x=v_usage_b.index, y=v_usage_b.values)
    plt.title("TRBV gene usage")
    plt.xticks(rotation=45)
    plt.show()
    print("TRBV usage:\n", v_usage_b)

if "TRBJ" in df.columns:
    j_usage_b = df["TRBJ"].value_counts()
    plt.figure(figsize=(10,5))
    sns.barplot(x=j_usage_b.index, y=j_usage_b.values)
    plt.title("TRBJ gene usage")
    plt.xticks(rotation=45)
    plt.show()
    print("TRBJ usage:\n", j_usage_b)

# Motifs β
df["motif_b"] = df["CDR3b"].str[-3:]
motif_counts_b = df["motif_b"].value_counts().head(10)
print("\nTop motifs β:\n", motif_counts_b)

# Diversity β
def shannon_diversity(sequences):
    counts = Counter(sequences)
    total = sum(counts.values())
    return -sum((c/total)*np.log(c/total) for c in counts.values())

print("Shannon diversity β:", shannon_diversity(df["CDR3b"]))

# =========================
# SECTION B: α-chain analysis
# =========================
if "CDR3a" in df.columns:
    df["cdr3a_length"] = df["CDR3a"].str.len()
    print("\nDistribution of CDR3α lengths:", df["cdr3a_length"])

    plt.figure(figsize=(8,5))
    sns.histplot(df["cdr3a_length"], bins=15, kde=True, color="blue")
    plt.title("CDR3α length distribution")
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.show()

if "TRAV" in df.columns:
    v_usage_a = df["TRAV"].value_counts()
    plt.figure(figsize=(10,5))
    sns.barplot(x=v_usage_a.index, y=v_usage_a.values)
    plt.title("TRAV gene usage")
    plt.xticks(rotation=45)
    plt.show()
    print("TRAV usage:\n", v_usage_a)

if "TRAJ" in df.columns:
    j_usage_a = df["TRAJ"].value_counts()
    plt.figure(figsize=(10,5))
    sns.barplot(x=j_usage_a.index, y=j_usage_a.values)
    plt.title("TRAJ gene usage")
    plt.xticks(rotation=45)
    plt.show()
    print("TRAJ usage:\n", j_usage_a)

# Motifs α
df["motif_a"] = df["CDR3a"].str[-3:]
motif_counts_a = df["motif_a"].value_counts().head(10)
print("\nTop motifs α:\n", motif_counts_a)

# Diversity α
print("Shannon diversity α:", shannon_diversity(df["CDR3a"]))

# =========================
# SECTION C: De-orphanization
# =========================
total = df.shape[0]
mapped = total   # كلهم mapped لأن الملف خاص بـEpitope واحد
orphans = total - mapped
rate = mapped / total * 100

print("\nDe-orphanization rate:")
print("Total:", total, "Mapped:", mapped, "Orphans:", orphans, "Rate: %.2f%%" % rate)

