import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# =========================
# SECTION 1: قراءة الخمس ملفات
# =========================
files = {
    "ELAGIGILTV":"VDJdb_ELAGIGILTV_TRA_TRB_Input.csv",
    "AAGIGILTV":"VDJdb_AAGIGILTV_TRA_TRB_Input.csv",
    "IMDQVPFSV":"VDJdb_IMDQVPFSV_TRA_TRB_Input.csv",
    "YLEPGPVTA":"VDJdb_YLEPGPVTA_TRA_TRB_Input.csv",
    "KTWGQYWQV":"VDJdb_KTWGQYWQV_TRA_TRB_Input.csv"
}

dfs = []
for epitope, fname in files.items():
    df = pd.read_csv(fname, sep="\t")
    df["epitope"] = epitope
    dfs.append(df)

all_df = pd.concat(dfs)

# =========================
# SECTION 2: جدول TCR–epitope
# =========================
tcr_map = all_df.groupby("CDR3b")["epitope"].unique()
print("\n=== TCR–Epitope mapping ===")
print(tcr_map)

# =========================
# SECTION 3: Cross-reactive TCRs
# =========================
cross_reactive = tcr_map[tcr_map.apply(len) > 1]
print("\n=== Cross-reactive TCRs (appear with >1 epitope) ===")
print(cross_reactive)

# =========================
# SECTION 4: Network graph
# =========================
G = nx.Graph()

# أضف nodes للـepitopes
for e in files.keys():
    G.add_node(e, bipartite="epitope", color="lightblue")

# أضف nodes للـTCRs وروابطهم
for tcr, eps in tcr_map.items():
    G.add_node(tcr, bipartite="tcr", color="lightgreen")
    for e in eps:
        G.add_edge(tcr, e)

# رسم الشبكة
colors = [G.nodes[n]["color"] for n in G.nodes()]
plt.figure(figsize=(14,10))
nx.draw(G, with_labels=False, node_color=colors, font_size=7, node_size=300)
plt.title("TCR–Epitope Network (VDJdb filtered data)")
plt.show()
