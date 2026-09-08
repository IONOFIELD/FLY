"""
Phase 1: sanity checks on the cached MaleCNS graph, plus locating
the neurons for the first cascade test (photoreceptors -> DNg13).

Run:  python inspect_graph.py
"""

import pandas as pd

neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")

# ------------------------------------------------ fix histamine sign
# Histamine is inhibitory in fly; modulatory amines get 0 for now.
sign_map = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "histamine": -1,
            "dopamine": 0, "octopamine": 0, "serotonin": 0}
edges["sign"] = edges["nt"].map(sign_map)
print("sign coverage:", edges["sign"].notna().mean().round(3))
edges.to_parquet("data/edges.parquet", index=False)

# ------------------------------------------------ basic graph stats
print("\nsuperclass counts:\n", neurons["superclass"].value_counts(dropna=False))
print("\nstatus counts:\n", neurons["status"].value_counts(dropna=False))

deg_out = edges.groupby("pre")["weight"].sum()
deg_in = edges.groupby("post")["weight"].sum()
print(f"\nneurons with any edge: {len(set(edges.pre) | set(edges.post)):,}")
print(f"median synapses out: {deg_out.median():.0f}, in: {deg_in.median():.0f}")
print(f"max synapses out: {deg_out.max():,}, in: {deg_in.max():,}")

# ------------------------------------------------ locate cascade endpoints
def find(pattern):
    m = neurons[neurons["type"].fillna("").str.contains(pattern, case=False, regex=True)]
    return m

photo = find(r"^R[1-8]")
print("\nphotoreceptor-like types:\n", photo["type"].value_counts().head(15))

dn = find(r"^DNg13")
print("\nDNg13 rows:\n", dn[["bodyId", "type", "instance", "somaSide", "predictedNt"]])

# ------------------------------------------------ shortest-path check
# Confirm a signed path exists from photoreceptors to DNg13 within a few hops.
import networkx as nx
G = nx.DiGraph()
G.add_weighted_edges_from(edges[["pre", "post", "weight"]].itertuples(index=False))
src = [b for b in photo["bodyId"] if b in G]
tgt = [b for b in dn["bodyId"] if b in G]
print(f"\nphotoreceptors in graph: {len(src):,}, DNg13 in graph: {len(tgt)}")
if src and tgt:
    lengths = []
    for t in tgt:
        try:
            lengths.append(nx.shortest_path_length(G, src[0], t))
        except nx.NetworkXNoPath:
            lengths.append(None)
    print("hops from first photoreceptor to each DNg13:", lengths)
