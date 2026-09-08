"""Shared helpers for benchmark scripts: type discovery, pulse protocol, readout scoring."""
import re
import numpy as np
import pandas as pd


def find_types(neurons, patterns):
    """Return list of distinct types matching any regex (case-insensitive), with counts."""
    t = neurons["type"].fillna("")
    hits = {}
    for pat in patterns:
        m = t[t.str.contains(pat, case=False, regex=True)]
        for k, v in m.value_counts().items():
            hits[k] = int(v)
    return hits


def present_types(neurons, wanted):
    have = set(neurons["type"].dropna().unique())
    return [w for w in wanted if w in have], [w for w in wanted if w not in have]


def pulse_protocol(model, rates_by_type, n_trials=10, pulse_ms=200, gap_ms=300, seed=0):
    """Drive stimulated types at given Hz for pulse_ms, n_trials times. Returns onsets."""
    rates = model.stim["type"].map(rates_by_type).fillna(0).values
    onsets = []
    for _ in range(n_trials):
        model.run(gap_ms)
        onsets.append(model.t_ms)
        model.set_stim_rates(rates)
        model.run(pulse_ms)
        model.set_stim_rates(np.zeros_like(rates))
    model.run(gap_ms)
    return onsets


def readout_rates(model, onsets, types, window_ms, meta):
    """Mean spikes per neuron of each type in the window after each onset, per side."""
    from .graph import infer_side
    side_all = infer_side(meta)
    sp = model.spike_frame()
    out = []
    for typ in types:
        ids = model.lif_ids[model.readout_index(typ)]
        if len(ids) == 0:
            continue
        side = side_all.loc[ids]
        for k, on in enumerate(onsets):
            w = sp[(sp.t_ms >= on) & (sp.t_ms < on + window_ms) & sp.bodyId.isin(ids)]
            for s in ["L", "R"]:
                n_s = (side == s).sum()
                if n_s:
                    out.append(dict(trial=k, type=typ, side=s,
                                    spikes_per_cell=w.bodyId.isin(side.index[side == s]).sum() / n_s,
                                    latency_ms=(w.t_ms.min() - on) if len(w) else np.nan))
    return pd.DataFrame(out)
