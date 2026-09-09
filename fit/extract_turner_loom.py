"""
Measured loom tuning from Turner, Krieger, Pang & Clandinin 2022 (eLife 11:e82587;
Dryad doi:10.5061/dryad.h44j0zpp8), read directly with h5py (no visanalysis).

Layout: datafiles/YYYY-MM-DD.hdf5 -> Flies/FlyN/epoch_runs/series_XXX/
   attrs: protocol_ID, num_epochs, pre_time, stim_time, tail_time, include_in_analysis
   aligned_response/glom/response  (gloms x frames), attrs mask_vals -> vpn_types.csv names
   acquisition attrs sample_period; stimulus_timing/frame_monitor (channels x 10 kHz samples)
   epochs/epoch_NNN attrs per trial (loom epochs carry 'Loom' in the radius trajectory)

Method: epoch onsets from the per-epoch `epoch_time` clock stamps relative to the series
`run_start_time` (verified evenly spaced: e.g. 100 epochs at 5.05 s, first 0.01 s, last
499.87 s). The frame monitor is a per-display-frame signal that is active almost continuously,
so rising-edge detection on it does not delimit epochs. Baseline = mean F in the pre window;
dF/F = (F - F0)/F0; response amplitude = peak dF/F in the stim window. Amplitudes are also
broken out by each epoch's `current_rv_ratio` (loom speed).

Behavioural state: each series carries a 50 Hz `binary_behavior` trace (the authors' own
walking classification, thresholded from frame-to-frame image RMSE) in its `behavior` group,
so no video download is needed. A loom trial is WALKING if the fly moves for >50% of the
stimulus window, STATIONARY if <10%, and is dropped otherwise; a series contributes to a
condition only with >= MIN_STATE_TRIALS trials in it. This matters: the two protocols agree on
glomerulus rank (Spearman 0.96) but differ in amplitude by up to 0.29 normalised units
(LC12 0.71 vs 1.00, LC17 1.00 vs 0.81), and the LoomingSpot series are the walking-fly days,
which is what Turner et al. 2022 is about.

Outputs results/turner2022/loom_tuning.json:
   amplitude[glom]  fly-averaged peak dF/F to loom, normalised to the max glom (all trials)
   amplitude_by_state[walking|stationary][glom]   same, split by behavioural state
   gain_sigma, gain_sigma_by_state                std of log(per-trial population gain)
   by_protocol      amplitudes computed separately for PanGlomSuite and LoomingSpot
Run:  python fit/extract_turner_loom.py [data/turner2022]
"""
import sys, json, glob, os, datetime
from pathlib import Path
import numpy as np
import pandas as pd
import h5py

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "data/turner2022")
OUT = Path("results/turner2022"); OUT.mkdir(parents=True, exist_ok=True)
vpn = pd.read_csv(ROOT / "template_brain" / "vpn_types.csv")
val2name = dict(zip(vpn.iloc[:, 0].astype(int), vpn["vpn_types"]))
INCLUDED = ["LC4", "LC6", "LC9", "LC11", "LC12", "LC15", "LC16", "LC17", "LC18", "LC21", "LC26", "LPLC1", "LPLC2"]
WALK_HI, WALK_LO, MIN_STATE_TRIALS = 0.5, 0.1, 5
LOOM_SET = ["LC6", "LC26", "LC16", "LPLC2", "LC4", "LPLC1", "LC9", "LC17", "LC12"]   # Turner Fig 5 loom-responsive


def _secs(t):
    h, m, rest = str(t).split(":")
    return int(h) * 3600 + int(m) * 60 + float(rest)


def epoch_table(series):
    """Per-epoch onset (s from series start), loom flag and r/v ratio, from the clock stamps."""
    t0 = _secs(series.attrs["run_start_time"])
    rows = []
    for k in sorted(series["epochs"].keys()):
        a = series["epochs"][k].attrs
        if "epoch_time" not in a:
            continue
        blob = " ".join(str(v) for v in a.values())
        rows.append(dict(onset=_secs(a["epoch_time"]) - t0,
                         is_loom=("Loom" in blob),
                         rv=float(a.get("current_rv_ratio", np.nan))))
    return pd.DataFrame(rows)


per_fly, per_trial_gain, sources, by_rv = [], [], [], []
by_state = {"walking": [], "stationary": []}
gain_state = {"walking": [], "stationary": []}
by_proto = {}


def walking_fraction(series, on, stim_end):
    """Fraction of the stimulus window the fly was moving, from the authors' binary trace."""
    b = series.get("behavior")
    if b is None or "binary_behavior" not in b or "frame_times" not in b:
        return np.nan
    ft = b["frame_times"][:]; bb = b["binary_behavior"][:]
    m = (ft[:len(bb)] >= on) & (ft[:len(bb)] < stim_end)
    return float(bb[m].mean()) if m.sum() >= 5 else np.nan
for path in sorted(glob.glob(str(ROOT / "datafiles" / "*.hdf5"))):
    with h5py.File(path, "r") as f:
        for fly in f["Flies"]:
            fa = f["Flies"][fly].attrs
            if str(fa.get("driver_1", "")) != "ChAT-T2A":
                continue
            for sname, s in f["Flies"][fly]["epoch_runs"].items():
                a = s.attrs
                proto = str(a.get("protocol_ID", ""))
                if proto not in ("LoomingSpot", "PanGlomSuite") or str(a.get("include_in_analysis", "True")) != "True":
                    continue
                if "aligned_response" not in s or "glom" not in s["aligned_response"]:
                    continue
                pre, stim, tail = float(a["pre_time"]), float(a["stim_time"]), float(a["tail_time"])
                n_exp = int(float(a["num_epochs"]))
                resp = s["aligned_response/glom/response"][:]
                names = [val2name.get(int(v), str(v)) for v in s["aligned_response/glom"].attrs["mask_vals"]]
                sp = float(s["acquisition"].attrs["sample_period"])
                fs = float(s["stimulus_timing"].attrs.get("sample_rate", 10000))
                fm = s["stimulus_timing/frame_monitor"][:]
                et = epoch_table(s)
                if not len(et):
                    print(f"  skip {os.path.basename(path)} {fly} {sname} ({proto}): no epoch timestamps")
                    continue
                if proto == "LoomingSpot":
                    et["is_loom"] = True
                t_frames = np.arange(resp.shape[1]) * sp
                peaks, peaks_rv = {}, {}
                state_peaks = {"walking": {}, "stationary": {}}
                state_gain = {"walking": [], "stationary": []}
                for on, lo, rv in zip(et.onset, et.is_loom, et.rv):
                    if not lo:
                        continue
                    wf = walking_fraction(s, on + pre, on + pre + stim)
                    state = ("walking" if wf > WALK_HI else "stationary" if wf < WALK_LO else None) \
                            if wf == wf else None
                    base = (t_frames >= on) & (t_frames < on + pre)
                    win = (t_frames >= on + pre) & (t_frames < on + pre + stim)
                    if base.sum() < 3 or win.sum() < 5:
                        continue
                    F0 = resp[:, base].mean(axis=1)
                    dff = (resp[:, win] - F0[:, None]) / np.where(F0 == 0, np.nan, F0)[:, None]
                    pk = np.nanmax(dff, axis=1)
                    for g, v in zip(names, pk):
                        peaks.setdefault(g, []).append(v)
                        if rv == rv and g in INCLUDED:
                            peaks_rv.setdefault((g, rv), []).append(v)
                        if state:
                            state_peaks[state].setdefault(g, []).append(v)
                if not peaks:
                    continue
                df = pd.DataFrame({g: pd.Series(v) for g, v in peaks.items()})
                df = df[[g for g in INCLUDED if g in df.columns]]
                key = f"{os.path.basename(path)}:{fly}:{sname}"
                per_fly.append(df.mean(axis=0).rename(key))
                by_proto.setdefault(proto, []).append(df.mean(axis=0).rename(key))
                for st, pk_map in state_peaks.items():
                    if not pk_map:
                        continue
                    sdf = pd.DataFrame({g: pd.Series(v) for g, v in pk_map.items()})
                    sdf = sdf[[g for g in INCLUDED if g in sdf.columns]]
                    if sdf.shape[0] < MIN_STATE_TRIALS:
                        continue
                    by_state[st].append(sdf.mean(axis=0).rename(key))
                    lc = [g for g in LOOM_SET if g in sdf.columns]
                    gnorm = sdf[lc] / sdf[lc].mean(axis=0).replace(0, np.nan)
                    gain_state[st].append(gnorm.mean(axis=1).values)
                loom_cols = [g for g in LOOM_SET if g in df.columns]
                norm = df[loom_cols] / df[loom_cols].mean(axis=0).replace(0, np.nan)
                per_trial_gain.append(norm.mean(axis=1).values)
                sources.append(dict(file=os.path.basename(path), fly=fly, series=sname, protocol=proto,
                                    loom_trials=int(df.shape[0]), gloms=int(df.shape[1]),
                                    rv_ratios=sorted({float(r) for r in et.rv[et.is_loom] if r == r})))
                for (g, rv), vals in peaks_rv.items():
                    by_rv.append(dict(glom=g, rv=rv, peak=float(np.nanmean(vals)), n=len(vals),
                                      series=f"{os.path.basename(path)}:{fly}:{sname}"))
                print(f"  {os.path.basename(path)} {fly} {sname} {proto:<13} {df.shape[0]:>3} loom trials, {df.shape[1]} gloms")

if not per_fly:
    raise SystemExit("no loom series extracted; check protocol names and photodiode detection")
amp = pd.concat(per_fly, axis=1)                       # gloms x series
mean_amp = amp.mean(axis=1); sem_amp = amp.std(axis=1) / np.sqrt(amp.shape[1])
normed = (mean_amp / mean_amp.max()).clip(lower=0)
gains = np.concatenate(per_trial_gain); gains = gains[np.isfinite(gains) & (gains > 0)]
out = dict(amplitude={k: round(float(v), 4) for k, v in normed.sort_values(ascending=False).items()},
           amplitude_sem_normalised={k: round(float(v / mean_amp.max()), 4) for k, v in sem_amp.items()},
           mean_peak_dff={k: round(float(v), 4) for k, v in mean_amp.items()},
           gain_sigma=round(float(np.std(np.log(gains))), 4), gain_n_trials=int(len(gains)),
           n_series=int(amp.shape[1]), n_flies=len({(s['file'], s['fly']) for s in sources}),
           behaviour_rule=dict(walking_gt=WALK_HI, stationary_lt=WALK_LO, min_trials=MIN_STATE_TRIALS,
                               source="per-series binary_behavior at 50 Hz (authors' classification)"),
           series=sources, method="photodiode epoch onsets; dF/F vs 1 s pre; peak in 3 s stim window",
           source="Turner, Krieger, Pang & Clandinin 2022, eLife 11:e82587; Dryad doi:10.5061/dryad.h44j0zpp8",
           extracted_on=datetime.date.today().isoformat())
# per-protocol amplitudes, reported because they differ
out["by_protocol"] = {}
for pr, cols in by_proto.items():
    a = pd.concat(cols, axis=1); mm = a.mean(axis=1)
    out["by_protocol"][pr] = dict(n_series=int(a.shape[1]),
                                  amplitude={k: round(float(v), 4) for k, v in (mm / mm.max()).items()})
# behavioural split
out["amplitude_by_state"], out["gain_sigma_by_state"], out["n_series_by_state"] = {}, {}, {}
for st, cols in by_state.items():
    if not cols:
        continue
    a = pd.concat(cols, axis=1); mm = a.mean(axis=1)
    out["amplitude_by_state"][st] = {k: round(float(v), 4) for k, v in (mm / mm.max()).items()}
    out["n_series_by_state"][st] = int(a.shape[1])
    g = np.concatenate(gain_state[st]) if gain_state[st] else np.array([])
    g = g[np.isfinite(g) & (g > 0)]
    out["gain_sigma_by_state"][st] = round(float(np.std(np.log(g))), 4) if len(g) > 20 else None
(OUT / "loom_tuning.json").write_text(json.dumps(out, indent=2))
if out["amplitude_by_state"]:
    st = pd.DataFrame(out["amplitude_by_state"]).dropna()
    print("\nnormalised amplitude by behavioural state:"); print(st.round(3).to_string())
    print("series per state:", out["n_series_by_state"], " gain sigma:", out["gain_sigma_by_state"])
if out.get("by_protocol"):
    pp = pd.DataFrame({k: v["amplitude"] for k, v in out["by_protocol"].items()}).dropna()
    print("\nnormalised amplitude by protocol:"); print(pp.round(3).to_string())
amp.to_csv(OUT / "loom_peak_dff_by_series.csv")
if by_rv:
    rvdf = pd.DataFrame(by_rv)
    rvdf.to_csv(OUT / "loom_peak_dff_by_rv.csv", index=False)
    piv = rvdf.pivot_table(index="glom", columns="rv", values="peak", aggfunc="mean")
    out["by_rv_ratio"] = {str(c): {k: round(float(v), 4) for k, v in piv[c].dropna().items()} for c in piv.columns}
    print("\npeak dF/F by loom r/v ratio (s):"); print(piv.round(3).to_string())
print("\nnormalised loom amplitude by glomerulus (fly-averaged):")
for k, v in out["amplitude"].items():
    print(f"  {k:<6} {v:.3f}")
print(f"trial gain sigma (log): {out['gain_sigma']}   from {out['gain_n_trials']} loom trials, {out['n_flies']} flies, {out['n_series']} series")
print(f"\nto activate:  cp {OUT/'loom_tuning.json'} data/loom_tuning_measured.json")
