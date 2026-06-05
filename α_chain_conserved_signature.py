import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ملفات لكل epitope
files = {
    "ELAGIGILTV": "VDJdb_ELAGIGILTV_TRA_TRB_Input.csv",
    "AAGIGILTV": "VDJdb_AAGIGILTV_TRA_TRB_Input.csv",
    "IMDQVPFSV": "VDJdb_IMDQVPFSV_TRA_TRB_Input.csv",
    "KTWGQYWQV": "VDJdb_KTWGQYWQV_TRA_TRB_Input.csv",
    "YLEPGPVTA": "VDJdb_YLEPGPVTA_TRA_TRB_Input.csv"
}

dfs = []
for ep, path in files.items():
    df = pd.read_csv(path, sep="\t")
    df["epitope"] = ep
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

# =========================
# TRAV usage across epitopes
trav_counts = all_df.groupby(["epitope","TRAV"]).size().unstack(fill_value=0)
trav_percent = trav_counts.div(trav_counts.sum(axis=1), axis=0) * 100

print("\n=== TRAV usage (%) per epitope ===")
print(trav_percent)

plt.figure(figsize=(12,6))
sns.heatmap(trav_percent, cmap="Blues", annot=True, fmt=".1f")
plt.title("TRAV usage across epitopes (%)")
plt.show()

# =========================
# TRAJ usage across epitopes
traj_counts = all_df.groupby(["epitope","TRAJ"]).size().unstack(fill_value=0)
traj_percent = traj_counts.div(traj_counts.sum(axis=1), axis=0) * 100

print("\n=== TRAJ usage (%) per epitope ===")
print(traj_percent)

plt.figure(figsize=(12,6))
sns.heatmap(traj_percent, cmap="Greens", annot=True, fmt=".1f")
plt.title("TRAJ usage across epitopes (%)")
plt.show()

# =========================
# Motif analysis (CDR3α suffix motifs)
all_df["motif_a"] = all_df["CDR3a"].str[-3:]
motif_counts = all_df.groupby(["epitope","motif_a"]).size().reset_index(name="count")

print("\n=== Top motifs in CDR3α per epitope ===")
for ep in all_df["epitope"].unique():
    subset = motif_counts[motif_counts["epitope"] == ep].sort_values("count", ascending=False).head(10)
    print(f"\n{ep}:\n", subset)

plt.figure(figsize=(12,6))
sns.barplot(data=motif_counts, x="motif_a", y="count", hue="epitope")
plt.title("CDR3α motifs across epitopes")
plt.xticks(rotation=45)
plt.show()
