"""
Classify taste-region sensory afferents by their signed influence on MN9.

For each sensory-superclass type T, influence on MN9 is
    sum over paths T -> M -> MN9 of  (w1*s1) * (w2*s2)   plus direct  w*s
normalised per neuron of T. Net-excitatory types are sugar-like candidates,
net-inhibitory types bitter-like (Gordon & Scott 2009: sugar excites, bitter
suppresses MN9). This is a wiring-based proxy for an annotation MaleCNS
lacks; it is reported, not asserted.

Run:  python benchmarks/diagnose_feeding.py
Writes results/feeding/afferent_influence.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph

DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")
edges["ws"] = edges["weight"] * edges["sign"]

mn9 = meta.index[meta["type"] == "MN9"]
into_mn9 = edges[edges.post.isin(mn9)]                       # M -> MN9
mid = into_mn9.groupby("pre").ws.sum()                         # signed weight of each M onto MN9
into_mid = edges[edges.post.isin(mid.index)]                   # T -> M
two_hop = into_mid.assign(w2=lambda d: d.post.map(mid)).eval("infl = ws * w2")
direct = into_mn9.assign(infl=lambda d: d.ws)

infl = pd.concat([two_hop[["pre", "infl"]], direct[["pre", "infl"]]])
infl["type"] = infl.pre.map(meta["type"])
infl["sc"] = infl.pre.map(meta["superclass"]).fillna("")
sens = infl[infl.sc.str.contains("sensory")]
tab = (sens.groupby("type").agg(influence=("infl", "sum"), n_pre=("pre", "nunique"))
       .assign(per_neuron=lambda d: d.influence / d.n_pre)
       .join(meta.groupby("type").size().rename("n_type")))
tab["nt"] = tab.index.map(lambda t: meta.loc[meta.type == t, "consensusNt"].mode().iat[0]
                          if (meta.type == t).any() else None)
tab = tab.sort_values("per_neuron", ascending=False)
out = Path("results/feeding"); out.mkdir(parents=True, exist_ok=True)
tab.to_csv(out / "afferent_influence.csv")

pd.set_option("display.width", 160)
print("\nNET EXCITATORY onto MN9 (sugar-like candidates):")
print(tab[tab.per_neuron > 0].head(12))
print("\nNET INHIBITORY onto MN9 (bitter-like candidates):")
print(tab[tab.per_neuron < 0].tail(12))
exc = list(tab[(tab.per_neuron > 0) & (tab.n_type >= 5)].head(6).index)
inh = list(tab[(tab.per_neuron < 0) & (tab.n_type >= 5)].tail(6).index)
print(f"\nsuggested drive:  python benchmarks/feeding.py 6 data {','.join(exc)} {','.join(inh)}")
