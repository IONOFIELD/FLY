"""
Why the model does not produce sugar-driven proboscis extension: a measured account.

All numbers below come from the MaleCNS graph and this model; the script writes
results/feeding/mechanism.json so the claim is reproducible and quotable.

M1 MN9's synaptic input is balanced: net signed weight near zero, so no uniform scaling
   of anything can tip it.
M2 Gustatory afferents (labellar bristle, taste peg, pharyngeal sensillum) contact MN9's
   inhibitory relays several times more strongly than its excitatory ones.
M3 Signed path products from those afferents to MN9 are negative at two hops and positive
   at three to five: the pathway is DISINHIBITORY.
M4 A silent network cannot express disinhibition, and a uniform tonic baseline that makes
   it non-silent removes specificity (the wind control drives MN9 as strongly as taste)
   and abolishes the escape response (benchmarks/fit_baseline.py).
Conclusion: this model class, with these transmitter predictions and no state-dependent
modulation, cannot produce the feeding response. Shiu et al. 2024 noted the same limitation
for inhibitory neurons in their 0 Hz-baseline model.

Run:  python benchmarks/feeding_mechanism.py [data_dir]
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns.graph import SIGN_MAP, GUSTATORY_MN9_DRIVING, GUSTATORY_SUBCLASSES

DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = Path("results/feeding"); OUT.mkdir(parents=True, exist_ok=True)
n = pd.read_parquet(f"{DATA}/neurons.parquet").set_index("bodyId")
e = pd.read_parquet(f"{DATA}/edges.parquet")
e["sgn"] = e["nt"].str.lower().map(SIGN_MAP).fillna(0)
mn9 = n.index[n.type == "MN9"]

# M1 input balance at MN9
inp = e[e.post.isin(mn9)]
exc_w = float(inp.loc[inp.sgn > 0, "weight"].sum()); inh_w = float(inp.loc[inp.sgn < 0, "weight"].sum())
M1 = dict(excitatory_synapses=exc_w, inhibitory_synapses=inh_w, net=exc_w - inh_w,
          net_fraction_of_total=round((exc_w - inh_w) / (exc_w + inh_w), 4))

# M2 where the afferents land
exc_src = set(inp.loc[inp.sgn > 0, "pre"]); inh_src = set(inp.loc[inp.sgn < 0, "pre"])
sub = n["subclass"].fillna("") if "subclass" in n else pd.Series("", index=n.index)
sets = {"screen_selected_subset": n.index[n.type.isin(GUSTATORY_MN9_DRIVING)],
        "all_gustatory_by_subclass": n.index[sub.isin(GUSTATORY_SUBCLASSES)]}
M2 = {}
for label, aff in sets.items():
    d = e[e.pre.isin(aff)]
    de = float(d[d.post.isin(exc_src)].weight.sum()); di = float(d[d.post.isin(inh_src)].weight.sum())
    M2[label] = dict(n_afferents=int(len(aff)), onto_MN9_excitatory=de, onto_MN9_inhibitory=di,
                     ratio=round(de / max(di, 1), 3))

# M3 signed path products, normalised by input weight
ei = e[e.sgn != 0]
ids = pd.Index(sorted(set(ei.pre) | set(ei.post))); idx = pd.Series(np.arange(len(ids)), index=ids)
W = csr_matrix((ei.weight.values * ei.sgn.values, (idx[ei.pre].values, idx[ei.post].values)),
               shape=(len(ids), len(ids)))
Wn = W.multiply(1.0 / np.maximum(np.abs(W).sum(axis=0), 1))
v = np.zeros(len(ids))
aff = [b for b in sets["screen_selected_subset"] if b in idx.index]
v[idx[aff].values] = 1.0
tgt = idx[[b for b in mn9 if b in idx.index]].values
M3 = {}
for hop in range(1, 6):
    v = np.asarray(Wn.T.dot(v)).ravel()
    M3[f"hop_{hop}"] = round(float(v[tgt].sum()), 5)

report = dict(
    M1_MN9_input_balance=M1,
    M2_afferent_targets=M2,
    M3_signed_path_products=M3,
    M4_uniform_manipulations_tested=[
        dict(manipulation="global synaptic scale (w_scale 0.25-1.0)", result="no value gives taste->MN9 with a quiet CNS",
             script="benchmarks/screen_afferents.py"),
        dict(manipulation="class-wise gain (sensory/relay/local, up to 3x)", result="MN9 silent at every combination",
             script="benchmarks/fit_gains.py"),
        dict(manipulation="regional gain on SEZ-intrinsic neurons (1.5-4x)",
             result="the SEZ population is net inhibitory onto the path under both a type-prefix and an "
                    "anatomical definition, so scaling amplifies inhibition with excitation",
             script="benchmarks/fit_regional.py"),
        dict(manipulation="uniform tonic depolarisation + membrane noise (3-6 mV)",
             result="MN9 responds but so does the wind control (specificity lost) and the escape response is abolished",
             script="benchmarks/fit_baseline.py")],
    conclusion=("taste -> MN9 in MaleCNS is disinhibitory; a silent network cannot express that, and no uniform "
                "manipulation restores it with specificity. Reproducing feeding needs cell-specific spontaneous "
                "activity, i.e. the state-dependent modulation this model omits (cf. Shiu et al. 2024, which "
                "reported the same limitation for inhibitory neurons at 0 Hz baseline)."),
    sign_map=SIGN_MAP)
(OUT / "mechanism.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
