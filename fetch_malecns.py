"""
Phase 0: pull the MaleCNS v1.0 graph from neuprint and cache it locally.

Run once:  python fetch_malecns.py
Outputs:   data/neurons.parquet   one row per neuron
           data/edges.parquet     one row per (pre, post) pair, weight >= MIN_WEIGHT

Requires:  export NEUPRINT_TOKEN="..."   (from neuprint.janelia.org account page)
"""

import os
import time
from pathlib import Path

import pandas as pd
from neuprint import Client, fetch_custom

# ---------------------------------------------------------------- settings
DATASET = "male-cns:v1.0"
MIN_WEIGHT = 5          # drop pre->post pairs with fewer than this many synapses
BATCH = 2000            # source neurons per edge query; lower if queries time out
OUT = Path("data")
OUT.mkdir(exist_ok=True)

token = os.environ.get("NEUPRINT_TOKEN")
if not token:
    raise SystemExit("Set NEUPRINT_TOKEN first: export NEUPRINT_TOKEN=\"...\"")

client = Client("neuprint.janelia.org", dataset=DATASET, token=token)
print("connected to", DATASET)

# ---------------------------------------------------------------- neurons
# Cypher returns null for properties a neuron lacks, so missing fields
# will not crash the query. We ask for every field we might want and
# inspect what actually comes back.
neuron_query = """
MATCH (n:Neuron)
RETURN n.bodyId        AS bodyId,
       n.type          AS type,
       n.instance      AS instance,
       n.status        AS status,
       n.predictedNt   AS predictedNt,
       n.consensusNt   AS consensusNt,
       n.somaLocation  AS somaLocation,
       n.somaSide      AS somaSide,
       n.superclass    AS superclass,
       n.pre           AS pre,
       n.post          AS post
"""
t0 = time.time()
neurons = fetch_custom(neuron_query, client=client)
print(f"neurons: {len(neurons):,} rows in {time.time()-t0:.0f}s")
print("columns present:", list(neurons.columns))
print("non-null counts:\n", neurons.notna().sum())

# somaLocation arrives as a dict like {'coordinates': [x, y, z]}; flatten it
def _xyz(v):
    if isinstance(v, dict) and "coordinates" in v:
        return v["coordinates"]
    return [None, None, None]

xyz = neurons["somaLocation"].apply(_xyz).tolist()
neurons[["x", "y", "z"]] = pd.DataFrame(xyz, index=neurons.index)
neurons = neurons.drop(columns=["somaLocation"])
neurons.to_parquet(OUT / "neurons.parquet", index=False)
print("wrote", OUT / "neurons.parquet")

# ---------------------------------------------------------------- edges
# One query for 125M synapses will time out, so batch by source bodyId.
edge_query = """
MATCH (a:Neuron)-[c:ConnectsTo]->(b:Neuron)
WHERE a.bodyId IN {ids} AND c.weight >= {minw}
RETURN a.bodyId AS pre, b.bodyId AS post, c.weight AS weight
"""
ids = neurons["bodyId"].tolist()
chunks = []
t0 = time.time()
for i in range(0, len(ids), BATCH):
    batch_ids = ids[i:i + BATCH]
    q = edge_query.format(ids=batch_ids, minw=MIN_WEIGHT)
    df = fetch_custom(q, client=client)
    chunks.append(df)
    done = min(i + BATCH, len(ids))
    print(f"  edges: {done:>7,}/{len(ids):,} sources, "
          f"{sum(len(c) for c in chunks):,} pairs, {time.time()-t0:.0f}s")

edges = pd.concat(chunks, ignore_index=True)

# ---------------------------------------------------------------- sign edges
# Fly rule: acetylcholine excitatory; GABA and glutamate inhibitory.
# Prefer consensusNt, fall back to predictedNt, else unknown.
nt = neurons.set_index("bodyId")
nt_col = nt["consensusNt"].fillna(nt["predictedNt"]).str.lower()
sign_map = {"acetylcholine": 1, "gaba": -1, "glutamate": -1}
edges["nt"] = edges["pre"].map(nt_col)
edges["sign"] = edges["nt"].map(sign_map)

n_total = len(edges)
n_unknown = edges["sign"].isna().sum()
print(f"\nedge pairs total:   {n_total:,}")
print(f"unknown NT sign:    {n_unknown:,} ({100*n_unknown/n_total:.1f}%)")
print("NT breakdown:\n", edges["nt"].value_counts(dropna=False))

edges.to_parquet(OUT / "edges.parquet", index=False)
print("wrote", OUT / "edges.parquet")
