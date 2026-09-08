"""
Stimulus protocols with provenance.

loom_protocol: population loom after Turner, Krieger, Pang & Clandinin 2022
(eLife 11:e82587). Loom-responsive glomeruli with a shared lognormal trial
gain. TUNING is an ordinal placeholder from their Fig 3A groups until the
authors' per-glomerulus loom dF/F replaces it.
"""
import numpy as np

LOOM_TUNING = {
    "LC4": 1.0, "LPLC2": 1.0, "LC6": 1.0, "LC16": 0.8,
    "LC26": 0.6, "LPLC1": 0.6, "LC9": 0.5, "LC17": 0.5, "LC12": 0.4,
}


def loom_protocol(model, n_trials=20, peak_hz=5.0, burst_ms=60, gap_ms=300,
                  gain_sigma=0.5, seed=0):
    """Run n_trials looms; return list of (onset_ms, gain)."""
    rng = np.random.default_rng(seed)
    tuning = model.stim["type"].map(LOOM_TUNING).fillna(0).values
    trials = []
    for _ in range(n_trials):
        model.run(gap_ms)
        gain = rng.lognormal(0, gain_sigma)
        trials.append((model.t_ms, gain))
        model.set_stim_rates(tuning * peak_hz * gain)
        model.run(burst_ms)
        model.set_stim_rates(np.zeros_like(tuning))
    model.run(gap_ms)
    return trials
