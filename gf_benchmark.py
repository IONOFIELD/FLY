"""
Tier 1 benchmark test: giant fiber escape circuit.

A loom is modelled as a 60 ms Poisson burst in LC4 and LPLC2.
We run N_TRIALS looms and report, per trial, GF and TTMn spike counts
and GF latency from burst onset. Then we repeat on a degree-preserving
rewired connectome; the response must disappear or it was not wiring.

Published physiology to compare against:
  GF spikes per loom: 1 to 2      GF latency: ~20 to 60 ms after loom onset
  TTMn follows GF within ~1 ms (electrical synapse, undercounted here)

Run:  python gf_benchmark.py
"""

import numpy as np
import pandas as pd
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor,
                    Network, ms, mV, Hz, second, defaultclock, prefs)

prefs.codegen.target = "cython"
defaultclock.dt = 0.1 * ms

STIM_TYPES = ["LC4", "LPLC2"]
READOUTS = ["DNp01", "TTMn"]
W_SCALE = 0.3               # from the sweep: stable, specific regime
BURST_HZ = 100
BURST_MS = 60
GAP_MS = 300                # silence between looms
N_TRIALS = 10
SEED = 0

V_REST, V_THRESH, V_RESET = -52 * mV, -45 * mV, -52 * mV
TAU_M, TAU_SYN, REFRACTORY = 20 * ms, 5 * ms, 2.2 * ms
W_BASE = 0.275 * mV

neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")
edges = edges[edges["sign"].notna() & (edges["sign"] != 0)].reset_index(drop=True)
meta = neurons.set_index("bodyId")


def rewire(df, seed):
    """Degree-preserving null: shuffle which postsynaptic neuron each edge lands on.
    Keeps every neuron's out-degree, total synapse counts and NT sign; destroys
    the specific wiring. Approximate (in-degree preserved in distribution only)."""
    rng = np.random.default_rng(seed)
    out = df.copy()
    out["post"] = rng.permutation(out["post"].values)
    out = out[out["pre"] != out["post"]]
    return out


def run(edge_df, label):
    stim_ids = sorted(set(neurons.loc[neurons["type"].isin(STIM_TYPES), "bodyId"]))
    stim_index = pd.Series(np.arange(len(stim_ids)), index=stim_ids)
    all_ids = pd.unique(pd.concat([edge_df["pre"], edge_df["post"]]))
    lif_ids = np.array([b for b in all_ids if b not in stim_index.index])
    lif_index = pd.Series(np.arange(len(lif_ids)), index=lif_ids)

    eqs = """
    dv/dt = (V_REST - v + g) / TAU_M : volt (unless refractory)
    dg/dt = -g / TAU_SYN : volt
    """
    G = NeuronGroup(len(lif_ids), eqs, threshold="v > V_THRESH", reset="v = V_RESET",
                    refractory=REFRACTORY, method="exact")
    G.v = V_REST

    rec = edge_df[edge_df["pre"].isin(lif_index.index) & edge_df["post"].isin(lif_index.index)]
    S_rec = Synapses(G, G, "w : volt", on_pre="g_post += w")
    S_rec.connect(i=lif_index[rec["pre"]].values, j=lif_index[rec["post"]].values)
    S_rec.w = (rec["weight"] * rec["sign"]).values * W_BASE * W_SCALE

    P = PoissonGroup(len(stim_ids), rates=0 * Hz)
    se = edge_df[edge_df["pre"].isin(stim_index.index) & edge_df["post"].isin(lif_index.index)]
    S_stim = Synapses(P, G, "w : volt", on_pre="g_post += w")
    S_stim.connect(i=stim_index[se["pre"]].values, j=lif_index[se["post"]].values)
    S_stim.w = (se["weight"] * se["sign"]).values * W_BASE * W_SCALE

    spikes = SpikeMonitor(G)
    net = Network(G, P, S_rec, S_stim, spikes)

    ro = {r: [lif_index[b] for b in meta.index[meta["type"] == r] if b in lif_index.index]
          for r in READOUTS}

    onsets = []
    for k in range(N_TRIALS):
        net.run(GAP_MS * ms)
        onsets.append(float(net.t / ms))
        P.rates = BURST_HZ * Hz
        net.run(BURST_MS * ms)
        P.rates = 0 * Hz
    net.run(GAP_MS * ms)

    t = np.array(spikes.t / ms)
    i = np.array(spikes.i)
    print(f"\n=== {label} ===  LIF {len(lif_ids):,}  synapses {len(rec):,}  "
          f"whole-CNS rate {len(i)/len(lif_ids)/float(net.t/second):.3f} Hz/neuron")
    print(f"{'trial':>5} {'GF spikes':>10} {'GF latency ms':>14} {'TTMn spikes':>12}")
    rows = []
    for k, on in enumerate(onsets):
        win = (t >= on) & (t < on + BURST_MS + 100)
        gf = np.isin(i[win], ro["DNp01"])
        tt = np.isin(i[win], ro["TTMn"])
        lat = (t[win][gf].min() - on) if gf.any() else np.nan
        rows.append(dict(trial=k, gf=gf.sum(), lat=lat, ttmn=tt.sum()))
        print(f"{k:>5} {gf.sum():>10} {lat:>14.1f} {tt.sum():>12}")
    df = pd.DataFrame(rows)
    print(f"mean GF spikes/loom {df.gf.mean():.1f}, mean latency {df.lat.mean():.1f} ms, "
          f"mean TTMn {df.ttmn.mean():.1f}")
    return df


real = run(edges, "REAL WIRING")
null = run(rewire(edges, SEED), "REWIRED NULL")
real.to_csv("data/gf_real.csv", index=False)
null.to_csv("data/gf_null.csv", index=False)
