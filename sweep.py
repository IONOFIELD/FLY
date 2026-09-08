"""
Phase 2b: find the regime where the cascade is stable AND specific.

For each (weight scale, stim Hz) we drive two stimulus sets and report:
  pop      whole-CNS rate, Hz per neuron
  frac     fraction of LIF neurons that fired at all
  GF/TTMn  mean rate of the two DNp01 (giant fiber) and two TTMn neurons
Specificity = GF rate under looming drive / GF rate under auditory drive.
Want: pop < ~1 Hz, frac < ~0.1, DNg13 visual >> DNg13 auditory.

Run:  python sweep.py
"""

import numpy as np
import pandas as pd
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor,
                    Network, ms, mV, Hz, second, defaultclock, prefs)

prefs.codegen.target = "cython"
defaultclock.dt = 0.1 * ms

STIM_SETS = {
    "visual":   ["LC4", "LPLC2"],
    "auditory": ["JO-A1", "JO-A2", "JO-B1_a", "JO-B2"],
}
READOUTS = ["DNp01", "TTMn"]
WEIGHT_SCALES = [1.0, 0.5, 0.25]        # multiplies 0.275 mV
STIM_RATES = [50, 100, 150]             # Hz
RUN_S = 0.5

V_REST, V_THRESH, V_RESET = -52 * mV, -45 * mV, -52 * mV
TAU_M, TAU_SYN, REFRACTORY = 20 * ms, 5 * ms, 2.2 * ms
W_BASE = 0.275 * mV

neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")
edges = edges[edges["sign"].notna() & (edges["sign"] != 0)]
meta = neurons.set_index("bodyId")

# all stimulated neurons from every set are Poisson sources; unused sets stay at 0 Hz
stim_ids = sorted(set(neurons.loc[neurons["type"].isin(sum(STIM_SETS.values(), [])), "bodyId"]))
stim_index = pd.Series(np.arange(len(stim_ids)), index=stim_ids)
all_ids = pd.unique(pd.concat([edges["pre"], edges["post"]]))
lif_ids = np.array([b for b in all_ids if b not in stim_index.index])
lif_index = pd.Series(np.arange(len(lif_ids)), index=lif_ids)

eqs = """
dv/dt = (V_REST - v + g) / TAU_M : volt (unless refractory)
dg/dt = -g / TAU_SYN : volt
"""
G = NeuronGroup(len(lif_ids), eqs, threshold="v > V_THRESH", reset="v = V_RESET",
                refractory=REFRACTORY, method="exact")
G.v = V_REST

rec = edges[edges["pre"].isin(lif_index.index) & edges["post"].isin(lif_index.index)]
S_rec = Synapses(G, G, "w : volt", on_pre="g_post += w")
S_rec.connect(i=lif_index[rec["pre"]].values, j=lif_index[rec["post"]].values)
rec_w = (rec["weight"] * rec["sign"]).values

P = PoissonGroup(len(stim_ids), rates=0 * Hz)
se = edges[edges["pre"].isin(stim_index.index) & edges["post"].isin(lif_index.index)]
S_stim = Synapses(P, G, "w : volt", on_pre="g_post += w")
S_stim.connect(i=stim_index[se["pre"]].values, j=lif_index[se["post"]].values)
stim_w = (se["weight"] * se["sign"]).values

spikes = SpikeMonitor(G)
net = Network(G, P, S_rec, S_stim, spikes)
net.store("fresh")
print(f"built: {len(lif_ids):,} LIF, {len(rec):,} recurrent, {len(se):,} stimulus synapses\n")

ro_idx = {r: [lif_index[b] for b in meta.index[meta["type"] == r] if b in lif_index.index] for r in READOUTS}
rates_per_set = {name: np.zeros(len(stim_ids)) for name in STIM_SETS}
for name, types in STIM_SETS.items():
    ids = neurons.loc[neurons["type"].isin(types), "bodyId"]
    rates_per_set[name][stim_index[ids].values] = 1.0

rows = []
for ws in WEIGHT_SCALES:
    for hz in STIM_RATES:
        res = {}
        for name in STIM_SETS:
            net.restore("fresh")
            S_rec.w = rec_w * W_BASE * ws
            S_stim.w = stim_w * W_BASE * ws
            P.rates = rates_per_set[name] * hz * Hz
            net.run(RUN_S * second)
            i = np.array(spikes.i)
            res[name] = dict(
                pop=len(i) / len(lif_ids) / RUN_S,
                frac=len(np.unique(i)) / len(lif_ids),
                dn=np.mean([(i == k).sum() / RUN_S for k in ro_idx["DNp01"]]),
                mn=np.mean([(i == k).sum() / RUN_S for k in ro_idx["TTMn"]]),
            )
        v, a = res["visual"], res["auditory"]
        spec = v["dn"] / a["dn"] if a["dn"] > 0 else float("inf")
        rows.append(dict(wscale=ws, hz=hz,
                         pop_vis=v["pop"], frac_vis=v["frac"], gf_vis=v["dn"], ttm_vis=v["mn"],
                         pop_aud=a["pop"], frac_aud=a["frac"], gf_aud=a["dn"], ttm_aud=a["mn"],
                         specificity=spec))
        print(f"w={ws:<5} hz={hz:<4} | vis pop={v['pop']:.3f} frac={v['frac']:.3f} DNg13={v['dn']:5.1f} "
              f"| aud pop={a['pop']:.3f} frac={a['frac']:.3f} DNg13={a['dn']:5.1f} | spec={spec:.1f}")

pd.DataFrame(rows).to_csv("data/sweep_results.csv", index=False)
print("\nwrote data/sweep_results.csv")
