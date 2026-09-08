"""
Giant fiber escape benchmark on MaleCNS v1.0, physiologically grounded version.

MODEL
  Leaky integrate-and-fire, current-based exponential synapses, parameters
  from Shiu et al. 2024 (Nature) with a global weight scale fitted on MaleCNS
  (see sweep.py; 0.275 mV/synapse is ~3x too strong when the VNC is included).
  Chemical synapse sign from MaleCNS consensus neurotransmitter prediction:
  ACh +, GABA -, Glu -, His - (fly), amines 0 (modulatory, not modelled).

ELECTRICAL SYNAPSES
  Invisible to EM connectomes. Modelled as two components on documented pairs:
    1. ohmic coupling   I = g_gap * (v_pre - v_post)          (subthreshold)
    2. spikelet         v_post += A_spikelet on pre spike      (AP is not in LIF)
  GF -> TTMn and GF -> PSI are one-to-one relays in vivo with ~0.8 ms latency
  (Tanouye & Wyman 1980; Allen, Godenschwege, Tanouye & Phelan 2006). A_spikelet
  is set so a single GF spike is suprathreshold; that is the physiological
  constraint, not a free parameter. Pairing is ipsilateral by somaSide.

LOOM STIMULUS
  Turner, Krieger, Pang & Clandinin 2022 (eLife 11:e82587) imaged 13 optic
  glomeruli simultaneously. Loom drives LC6, LC26, LC16, LPLC2, LC4, LPLC1,
  LC9, LC17, LC12 with type-specific amplitude, a large shared trial-to-trial
  gain factor (their Fig 4), and no walking suppression on loom channels
  (their Fig 5E,F). We drive those types as Poisson sources with rates
  = PEAK_HZ * tuning[type] * gain_trial, gain_trial ~ lognormal.
  PLACEHOLDER: tuning[] below is a coarse ordinal from their Fig 3A groups.
  Replace with mean loom dF/F per glomerulus from the authors.
  PLACEHOLDER: dF/F -> spikes/s mapping unknown; PEAK_HZ anchored so GF fires
  1 to 2 spikes per loom (von Reyn et al. 2014; Ache et al. 2019).

PASS CRITERIA (published physiology)
  GF: 1 to 2 spikes per loom, hit rate near 1.   TTMn: 1 spike per GF spike,
  lag < 1.5 ms.   Whole-CNS rate stays < 0.1 Hz/neuron.   Rewired null: ~0.

Run:  python gf_physio.py
"""

import numpy as np
import pandas as pd
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor,
                    Network, ms, mV, Hz, second, defaultclock, prefs, seed)

prefs.codegen.target = "cython"
defaultclock.dt = 0.1 * ms
SEED = 0
seed(SEED)
rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------- parameters
W_SCALE = 0.3
V_REST, V_THRESH, V_RESET = -52 * mV, -45 * mV, -52 * mV
TAU_M, TAU_SYN, REFRACTORY = 20 * ms, 5 * ms, 2.2 * ms
W_BASE = 0.275 * mV

G_GAP = 0.2                 # dimensionless ohmic coupling (fraction of leak)
A_SPIKELET = 9 * mV         # > 7 mV rest-to-threshold gap: one-to-one relay

PEAK_HZ = 5                 # rate of the most loom-sensitive type at gain 1
BURST_MS = 60
GAP_MS = 300
N_TRIALS = 20
GAIN_SIGMA = 0.5            # lognormal sigma of shared trial gain (Turner Fig 4)

# PLACEHOLDER ordinal tuning from Turner et al. 2022 Fig 3A; replace with data
LOOM_TUNING = {
    "LC4": 1.0, "LPLC2": 1.0, "LC6": 1.0, "LC16": 0.8,
    "LC26": 0.6, "LPLC1": 0.6, "LC9": 0.5, "LC17": 0.5, "LC12": 0.4,
}

ELECTRICAL = [   # (pre type, post type, source)
    ("DNp01", "TTMn", "Allen et al. 2006; Tanouye & Wyman 1980"),
    ("DNp01", "PSI",  "Allen et al. 2006; Phelan et al. 2008"),
]
READOUTS = ["DNp01", "TTMn", "PSI"]

# ---------------------------------------------------------------- data
neurons = pd.read_parquet("data/neurons.parquet")
edges = pd.read_parquet("data/edges.parquet")
edges = edges[edges["sign"].notna() & (edges["sign"] != 0)].reset_index(drop=True)
meta = neurons.set_index("bodyId")


def electrical_pairs():
    """Ipsilateral pre->post pairs for each documented electrical synapse."""
    pairs = []
    for pre_t, post_t, src in ELECTRICAL:
        pres = meta[meta["type"] == pre_t]
        posts = meta[meta["type"] == post_t]
        for b, row in pres.iterrows():
            same = posts[posts["somaSide"] == row["somaSide"]]
            q = same.index[0] if len(same) else posts.index[0]
            pairs.append((b, q))
            print(f"electrical {pre_t} {b} ({row['somaSide']}) -> {post_t} {q}  [{src}]")
    return pairs


def rewire(df):
    out = df.copy()
    out["post"] = rng.permutation(out["post"].values)
    return out[out["pre"] != out["post"]]


def run(edge_df, label, with_electrical=True):
    stim_types = list(LOOM_TUNING)
    stim = neurons[neurons["type"].isin(stim_types)][["bodyId", "type"]]
    stim_index = pd.Series(np.arange(len(stim)), index=stim["bodyId"].values)
    tuning_vec = stim["type"].map(LOOM_TUNING).values

    all_ids = pd.unique(pd.concat([edge_df["pre"], edge_df["post"]]))
    lif_ids = np.array([b for b in all_ids if b not in stim_index.index])
    lif_index = pd.Series(np.arange(len(lif_ids)), index=lif_ids)

    eqs = """
    dv/dt = (V_REST - v + g + I_gap) / TAU_M : volt (unless refractory)
    dg/dt = -g / TAU_SYN : volt
    I_gap : volt
    """
    G = NeuronGroup(len(lif_ids), eqs, threshold="v > V_THRESH", reset="v = V_RESET",
                    refractory=REFRACTORY, method="euler")
    G.v = V_REST

    rec = edge_df[edge_df["pre"].isin(lif_index.index) & edge_df["post"].isin(lif_index.index)]
    S_rec = Synapses(G, G, "w : volt", on_pre="g_post += w")
    S_rec.connect(i=lif_index[rec["pre"]].values, j=lif_index[rec["post"]].values)
    S_rec.w = (rec["weight"] * rec["sign"]).values * W_BASE * W_SCALE

    P = PoissonGroup(len(stim), rates=0 * Hz)
    se = edge_df[edge_df["pre"].isin(stim_index.index) & edge_df["post"].isin(lif_index.index)]
    S_stim = Synapses(P, G, "w : volt", on_pre="g_post += w")
    S_stim.connect(i=stim_index[se["pre"]].values, j=lif_index[se["post"]].values)
    S_stim.w = (se["weight"] * se["sign"]).values * W_BASE * W_SCALE

    objs = [G, P, S_rec, S_stim]
    if with_electrical:
        pairs = [(lif_index[a], lif_index[b]) for a, b in electrical_pairs()
                 if a in lif_index.index and b in lif_index.index]
        S_gap = Synapses(G, G,
                         "I_gap_post = G_GAP * (v_pre - v_post) : volt (summed)",
                         on_pre="v_post += A_SPIKELET")
        S_gap.connect(i=[a for a, _ in pairs], j=[b for _, b in pairs])
        objs.append(S_gap)

    spikes = SpikeMonitor(G)
    objs.append(spikes)
    net = Network(*objs)

    ro = {r: [lif_index[b] for b in meta.index[meta["type"] == r] if b in lif_index.index]
          for r in READOUTS}

    onsets, gains = [], []
    for k in range(N_TRIALS):
        net.run(GAP_MS * ms)
        gain = rng.lognormal(0, GAIN_SIGMA)
        gains.append(gain)
        onsets.append(float(net.t / ms))
        P.rates = tuning_vec * PEAK_HZ * gain * Hz
        net.run(BURST_MS * ms)
        P.rates = 0 * Hz
    net.run(GAP_MS * ms)

    t = np.array(spikes.t / ms)
    i = np.array(spikes.i)
    print(f"\n=== {label} ===  LIF {len(lif_ids):,}  chem synapses {len(rec):,}  "
          f"stim neurons {len(stim):,}  whole-CNS {len(i)/len(lif_ids)/float(net.t/second):.4f} Hz/neuron")
    print(f"{'trial':>5} {'gain':>6} {'GF/cell':>8} {'GF lat':>7} {'TTMn/cell':>10} {'TTMn lag':>9} {'PSI/cell':>9}")
    rows = []
    for k, on in enumerate(onsets):
        win = (t >= on) & (t < on + BURST_MS + 100)
        tw, iw = t[win], i[win]
        gf, tt, ps = (np.isin(iw, ro[r]) for r in ["DNp01", "TTMn", "PSI"])
        lat = tw[gf].min() - on if gf.any() else np.nan
        tlag = tw[tt].min() - tw[gf].min() if (gf.any() and tt.any()) else np.nan
        r = dict(trial=k, gain=gains[k], gf=gf.sum() / len(ro["DNp01"]), lat=lat,
                 ttmn=tt.sum() / len(ro["TTMn"]), tlag=tlag, psi=ps.sum() / max(len(ro["PSI"]), 1))
        rows.append(r)
        print(f"{k:>5} {r['gain']:>6.2f} {r['gf']:>8.1f} {lat:>7.1f} {r['ttmn']:>10.1f} {tlag:>9.2f} {r['psi']:>9.1f}")
    df = pd.DataFrame(rows)
    print(f"GF {df.gf.mean():.2f}/cell (hit {(df.gf>0).mean():.2f}), latency {df.lat.mean():.1f} ms | "
          f"TTMn {df.ttmn.mean():.2f}/cell, lag {df.tlag.mean():.2f} ms | PSI {df.psi.mean():.2f}/cell")
    return df


real = run(edges, "REAL WIRING + electrical")
chem = run(edges, "REAL WIRING, chemical only", with_electrical=False)
null = run(rewire(edges), "REWIRED NULL + electrical")
for name, df in [("real", real), ("chem_only", chem), ("null", null)]:
    df.to_csv(f"data/gf_physio_{name}.csv", index=False)
print("\nwrote data/gf_physio_*.csv")
