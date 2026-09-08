"""
Export simulated sessions in a brainsets-style HDF5 layout for POYO / POYO+ (torch_brain).

Each session: one stimulus family (loom / sound / taste), N trials with per-trial condition
labels, spikes from a recorded subset of neurons, and a units table carrying connectome
features that could seed unit embeddings:
    type, superclass, subclass, side, in_synapses, out_synapses, in_degree, out_degree,
    synaptic_depth_from_stimulus (BFS hops), and a stable unit_id (MaleCNS bodyId).

Layout (per file, groups/datasets):
    /spikes/timestamps  float64 [s]      /spikes/unit_index int32
    /units/id, /units/<feature>...      /trials/start, /trials/end, /trials/condition
    attrs: session_id, subject_id="malecns_v1.0_lif", brainset="flycns_sim", sampling=event
This follows the structure brainsets uses (temporaldata Data: spikes, units, trials,
domain); verify field names against the installed brainsets version before training.

Run:  python export_brainsets.py [n_trials] [out_dir]
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from flycns import load_graph, CNSModel, LIFParams
from flycns.protocols import LOOM_TUNING
from flycns.bench import find_types, pulse_protocol
from flycns.graph import GUSTATORY_MN9_DRIVING

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "results/brainsets"); OUT.mkdir(parents=True, exist_ok=True)
neurons, edges = load_graph("data")
meta = neurons.set_index("bodyId")
jo = [t for t in find_types(neurons, [r"^JO-A", r"^JO-B"]) if not t.endswith("unclear")]
taste = [t for t in GUSTATORY_MN9_DRIVING if (neurons.type == t).any()]

SESSIONS = {
    "loom":  dict(stim=list(LOOM_TUNING), conditions={"weak": 0.5, "medium": 1.0, "strong": 2.0}),
    "sound": dict(stim=jo,               conditions={"soft": 50, "loud": 150}),
    "taste": dict(stim=taste,            conditions={"low": 50, "high": 150}),
}
# recorded population: everything within 3 synapses downstream of any stimulus type (capped)
import networkx as nx
G = nx.DiGraph(); G.add_edges_from(edges[["pre", "post"]].itertuples(index=False))
deg_in = edges.groupby("post").weight.sum(); deg_out = edges.groupby("pre").weight.sum()
nin = edges.groupby("post").size(); nout = edges.groupby("pre").size()

for name, cfg in SESSIONS.items():
    stim_ids = set(neurons.loc[neurons.type.isin(cfg["stim"]), "bodyId"])
    depth = {}
    frontier = {b: 0 for b in stim_ids if b in G}
    seen = dict(frontier)
    for hop in range(1, 4):
        nxt = {}
        for b in frontier:
            for s in G.successors(b):
                if s not in seen:
                    seen[s] = hop; nxt[s] = hop
        frontier = nxt
    recorded = [b for b, h in seen.items() if h >= 1]
    rng = np.random.default_rng(0)
    if len(recorded) > 3000:
        recorded = list(rng.choice(recorded, 3000, replace=False))
    recorded = sorted(recorded)
    print(f"[{name}] {len(stim_ids)} stimulated, recording {len(recorded)} downstream neurons")

    m = CNSModel(neurons, edges, cfg["stim"], LIFParams(), electrical=True)
    tuning = m.stim["type"].map(LOOM_TUNING).fillna(0).values if name == "loom" else None
    trials = []
    for k in range(N):
        cond = list(cfg["conditions"])[k % len(cfg["conditions"])]
        val = cfg["conditions"][cond]
        m.run(300); t0 = m.t_ms
        rates = tuning * 5.0 * val if name == "loom" else np.full(len(m.stim), float(val))
        m.set_stim_rates(rates); m.run(60 if name == "loom" else 200)
        m.set_stim_rates(np.zeros(len(m.stim)))
        trials.append(dict(start=t0 / 1000, end=(m.t_ms + 200) / 1000, condition=cond, value=val))
    m.run(300)
    sp = m.spike_frame(); sp = sp[sp.bodyId.isin(recorded)]
    uidx = pd.Series(np.arange(len(recorded)), index=recorded)

    f = OUT / f"flycns_sim_{name}.h5"
    with h5py.File(f, "w") as h:
        h.attrs.update(dict(session_id=f"flycns_sim_{name}", subject_id="malecns_v1.0_lif", brainset="flycns_sim",
                            model_params=json.dumps(LIFParams().__dict__), stimulus_types=json.dumps(cfg["stim"]),
                            note="simulated; connectome-constrained LIF; see flycns RESULTS.md"))
        g = h.create_group("spikes")
        g.create_dataset("timestamps", data=(sp.t_ms.values / 1000.0).astype("f8"))
        g.create_dataset("unit_index", data=uidx.loc[sp.bodyId].values.astype("i4"))
        u = h.create_group("units")
        u.create_dataset("id", data=np.array([f"malecns:{b}" for b in recorded], dtype="S"))
        u.create_dataset("body_id", data=np.array(recorded, dtype="i8"))
        for col in ["type", "superclass", "subclass", "somaSide", "consensusNt"]:
            if col in meta:
                u.create_dataset(col, data=np.array(meta.loc[recorded, col].fillna("").astype(str), dtype="S"))
        u.create_dataset("in_synapses", data=deg_in.reindex(recorded).fillna(0).values.astype("f4"))
        u.create_dataset("out_synapses", data=deg_out.reindex(recorded).fillna(0).values.astype("f4"))
        u.create_dataset("in_degree", data=nin.reindex(recorded).fillna(0).values.astype("f4"))
        u.create_dataset("out_degree", data=nout.reindex(recorded).fillna(0).values.astype("f4"))
        u.create_dataset("synaptic_depth_from_stimulus", data=np.array([seen[b] for b in recorded], dtype="i2"))
        for c in ["x", "y", "z"]:
            u.create_dataset(f"soma_{c}", data=meta.loc[recorded, c].fillna(np.nan).values.astype("f4"))
        tr = h.create_group("trials")
        for key, dt in [("start", "f8"), ("end", "f8"), ("value", "f4")]:
            tr.create_dataset(key, data=np.array([t[key] for t in trials], dtype=dt))
        tr.create_dataset("condition", data=np.array([t["condition"] for t in trials], dtype="S"))
        h.create_dataset("domain/start", data=np.array([0.0])); h.create_dataset("domain/end", data=np.array([m.t_ms / 1000.0]))
    print(f"  wrote {f}: {len(sp):,} spikes, {len(recorded)} units, {N} trials, {m.t_ms/1000:.1f} s")
print("\nnext: pip install brainsets torch_brain; load with h5py or adapt brainsets' Data.from_hdf5 and verify field names")
