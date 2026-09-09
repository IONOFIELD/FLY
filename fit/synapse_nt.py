"""
Per-synapse neurotransmitter signing and a confidence-weighted edge table.

The model currently signs every synapse of a neuron with that neuron's aggregate
consensus transmitter. MaleCNS also publishes per-T-bar prediction probabilities
(tbar-neurotransmitters-male-cns-v1.0.feather), which let us:

  1. weight each connection by the classifier's confidence, so a 0.51 call does not
     count the same as a 0.99 one  (edges gain a `conf` column; the model can scale
     weights by it, or leave it unused);
  2. flag neurons whose T-bars disagree with their aggregate label, which is a
     data-driven replacement for the arbitrary 5%/10% random sign-flip perturbation.

Writes data/edges_nt.parquet (edges + mean per-synapse confidence + fraction of
T-bars agreeing with the neuron's consensus) and results/nt_confidence.json.

Run:  python fit/synapse_nt.py           (needs data/bulk/tbar-neurotransmitters.feather)
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns.graph import SIGN_MAP

TB = Path("data/bulk/tbar-neurotransmitters.feather")
if not TB.exists():
    raise SystemExit("missing data/bulk/tbar-neurotransmitters.feather; run: python fetch_bulk.py synapse-nt")
neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")
meta = neurons.set_index("bodyId")

print("reading per-T-bar predictions (2.7 GB, a minute or two)...")
tb = pd.read_feather(TB)
print("columns:", list(tb.columns)[:12])
# columns are nt_<transmitter>_prob in MaleCNS v1.0
nt_cols = [c for c in tb.columns if c.startswith("nt_") and c.endswith("_prob")]
nt_names = [c[3:-5].lower() for c in nt_cols]
body_col = next((c for c in tb.columns if c.lower() in ("body", "bodyid", "body_pre")), None)
print(f"transmitters: {nt_names}")
if not nt_cols or body_col is None:
    raise SystemExit(f"unexpected columns; found nt={nt_cols} body={body_col}. Inspect and adjust.")

p = tb[nt_cols].to_numpy()
best = np.argmax(p, axis=1)
conf = p[np.arange(len(p)), best]
tb_small = pd.DataFrame({"body": tb[body_col].values, "pred": [nt_names[i] for i in best], "conf": conf})

# The `body` column covers every EM segment (~1.8M), most of them unproofread fragments
# with no consensus label; restrict to the annotated neurons that carry the graph.
n_all = len(tb_small := tb_small if "tb_small" in dir() else None) if False else None
cons = meta["consensusNt"].str.lower()
n_before = tb_small["body"].nunique()
tb_small = tb_small[tb_small["body"].isin(cons.dropna().index)]
print(f"restricting to annotated neurons: {tb_small['body'].nunique():,} of {n_before:,} bodies")
tb_small["consensus"] = tb_small["body"].map(cons)
tb_small["agree"] = tb_small["pred"].values == tb_small["consensus"].fillna("").values
agg = tb_small.groupby("body").agg(mean_conf=("conf", "mean"), n_tbars=("conf", "size"),
                                   agree=("agree", "mean"))
edges2 = edges.copy()
edges2["nt_conf"] = edges2["pre"].map(agg["mean_conf"])
edges2["nt_agree"] = edges2["pre"].map(agg["agree"])
edges2.to_parquet("data/edges_nt.parquet", index=False)

# how much of the graph rests on low-confidence calls
out = {"n_neurons_with_tbars": int(len(agg)),
       "mean_confidence": round(float(agg.mean_conf.mean()), 4),
       "frac_neurons_below_0.6": round(float((agg.mean_conf < 0.6).mean()), 4),
       "frac_neurons_tbars_disagree_with_consensus": round(float((agg.agree < 0.5).mean()), 4),
       "edge_weight_fraction_from_low_conf": round(float(
           edges2.loc[edges2.nt_conf < 0.6, "weight"].sum() / edges2.weight.sum()), 4),
       "edge_weight_fraction_from_disagreeing_neurons": round(float(
           edges2.loc[edges2.nt_agree < 0.5, "weight"].sum() / edges2.weight.sum()), 4),
       "transmitters": nt_names,
       "note_bodies": ("statistics are over annotated neurons only; the raw file covers ~1.8M EM "
                       "segments, mostly unproofread fragments with no consensus label"),
       "note": ("per-T-bar predictions vs the aggregate consensus label used for signing; "
                "use nt_agree to target sign perturbations at genuinely uncertain neurons "
                "instead of a random fraction")}
Path("results").mkdir(exist_ok=True)
Path("results/nt_confidence.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
print("\nwrote data/edges_nt.parquet (edges + nt_conf + nt_agree)")
