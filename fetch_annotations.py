"""
Add sensory annotations to data/neurons.parquet without refetching edges:
subclass (e.g. 'pharyngeal sensillum', 'labellar bristle', 'taste peg') and
entryNerve (e.g. 'MxLbN', 'aPhN', 'AN'). Safe to rerun.

Run:  python fetch_annotations.py
"""
import os
import pandas as pd
from neuprint import Client, fetch_custom

token = os.environ.get("NEUPRINT_TOKEN")
if not token:
    raise SystemExit('Set NEUPRINT_TOKEN first: export NEUPRINT_TOKEN="..."')
c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=token)
q = """MATCH (n:Neuron) RETURN n.bodyId AS bodyId, n.subclass AS subclass,
       n.entryNerve AS entryNerve, n.exitNerve AS exitNerve"""
ann = fetch_custom(q, client=c)
n = pd.read_parquet("data/neurons.parquet")
for col in ["subclass", "entryNerve", "exitNerve"]:
    if col in n:
        n = n.drop(columns=col)
n = n.merge(ann, on="bodyId", how="left")
n.to_parquet("data/neurons.parquet", index=False)
print("annotated:", n.subclass.notna().sum(), "with subclass;", n.entryNerve.notna().sum(), "with entryNerve")
print(n.subclass.value_counts().head(20))
