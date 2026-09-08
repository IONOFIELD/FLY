"""
Phase 2: leaky integrate-and-fire on the full MaleCNS graph (Brian2).

Parameters follow Shiu et al. Nature 2024 (FlyWire LIF model).
Stimulated neurons are replaced by Poisson spike sources at STIM_HZ.

Run:  python simulate.py
Edit STIM_TYPES and READOUT_TYPES to change the experiment.
"""

import numpy as np
import pandas as pd
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor,
                    Network, ms, mV, Hz, second, defaultclock, prefs)

prefs.codegen.target = "cython"     # falls back to numpy if no compiler
defaultclock.dt = 0.1 * ms

# ------------------------------------------------ experiment settings
STIM_TYPES = ["JO-A1", "JO-A2", "JO-B1_a", "JO-B2"]              # types driven as Poisson sources
STIM_HZ = 150                        # firing rate of stimulated neurons
READOUT_TYPES = ["DNg13"]           # types whose spikes we report
BASELINE_S = 0.5                    # seconds with stimulus off
STIM_S = 0.5                        # seconds with stimulus on

# ------------------------------------------------ Shiu et al. parameters
V_REST = -52 * mV
V_THRESH = -45 * mV
V_RESET = -52 * mV
TAU_M = 20 * ms
TAU_SYN = 5 * ms
REFRACTORY = 2.2 * ms
W_PER_SYNAPSE = 0.275 * mV

# ------------------------------------------------ load graph
neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")
edges = edges[edges["sign"].notna() & (edges["sign"] != 0)]   # drop unknown and modulatory

stim_ids = set(neurons.loc[neurons["type"].isin(STIM_TYPES), "bodyId"])
print(f"stimulating {len(stim_ids):,} neurons of types {STIM_TYPES}")

# everything with an edge that is not a stimulated neuron becomes an LIF unit
all_ids = pd.unique(pd.concat([edges["pre"], edges["post"]]))
lif_ids = np.array([b for b in all_ids if b not in stim_ids])
lif_index = pd.Series(np.arange(len(lif_ids)), index=lif_ids)
stim_ids = np.array(sorted(stim_ids))
stim_index = pd.Series(np.arange(len(stim_ids)), index=stim_ids)
print(f"LIF neurons: {len(lif_ids):,}")

# ------------------------------------------------ build model
eqs = """
dv/dt = (V_REST - v + g) / TAU_M : volt (unless refractory)
dg/dt = -g / TAU_SYN : volt
"""
G = NeuronGroup(len(lif_ids), eqs, threshold="v > V_THRESH", reset="v = V_RESET",
                refractory=REFRACTORY, method="exact")
G.v = V_REST

# recurrent synapses among LIF neurons
rec = edges[edges["pre"].isin(lif_index.index) & edges["post"].isin(lif_index.index)]
S_rec = Synapses(G, G, "w : volt", on_pre="g_post += w")
S_rec.connect(i=lif_index[rec["pre"]].values, j=lif_index[rec["post"]].values)
S_rec.w = (rec["weight"] * rec["sign"]).values * W_PER_SYNAPSE
print(f"recurrent synapses: {len(rec):,}")

# stimulus sources, silent during baseline then Poisson at STIM_HZ
P = PoissonGroup(len(stim_ids), rates=0 * Hz)
stim_edges = edges[edges["pre"].isin(stim_index.index) & edges["post"].isin(lif_index.index)]
S_stim = Synapses(P, G, "w : volt", on_pre="g_post += w")
S_stim.connect(i=stim_index[stim_edges["pre"]].values, j=lif_index[stim_edges["post"]].values)
S_stim.w = (stim_edges["weight"] * stim_edges["sign"]).values * W_PER_SYNAPSE
print(f"stimulus synapses: {len(stim_edges):,}")

spikes = SpikeMonitor(G)
net = Network(G, P, S_rec, S_stim, spikes)

# ------------------------------------------------ run
print("running baseline...")
net.run(BASELINE_S * second, report="text")
P.rates = STIM_HZ * Hz
print("running stimulus...")
net.run(STIM_S * second, report="text")

# ------------------------------------------------ readout
t = np.array(spikes.t / second)
i = np.array(spikes.i)
base = t < BASELINE_S
stim = ~base
print(f"\npopulation rate baseline: {base.sum()/len(lif_ids)/BASELINE_S:.3f} Hz per neuron")
print(f"population rate stimulus: {stim.sum()/len(lif_ids)/STIM_S:.3f} Hz per neuron")

meta = neurons.set_index("bodyId")
for typ in READOUT_TYPES:
    for b in meta.index[meta["type"] == typ]:
        if b not in lif_index.index:
            continue
        k = lif_index[b]
        nb = ((i == k) & base).sum() / BASELINE_S
        ns = ((i == k) & stim).sum() / STIM_S
        print(f"{typ} {b} ({meta.at[b, 'somaSide']}): baseline {nb:.1f} Hz, stimulus {ns:.1f} Hz")

# types most changed by the stimulus
counts = pd.DataFrame({"idx": i, "stim": stim})
counts["bodyId"] = lif_ids[counts["idx"]]
counts["type"] = counts["bodyId"].map(meta["type"])
delta = (counts[counts.stim].groupby("type").size() / STIM_S
         - counts[~counts.stim].groupby("type").size().reindex(
             counts[counts.stim].groupby("type").size().index, fill_value=0) / BASELINE_S)
print("\ntop 20 types by added spikes per second during stimulus:")
print(delta.sort_values(ascending=False).head(20))
