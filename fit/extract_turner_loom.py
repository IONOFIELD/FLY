"""
Measured loom tuning from Turner, Krieger, Pang & Clandinin 2022 (eLife 11:e82587;
Dryad doi:10.5061/dryad.h44j0zpp8), read directly with h5py (no visanalysis).

Layout: datafiles/YYYY-MM-DD.hdf5 -> Flies/FlyN/epoch_runs/series_XXX/
   attrs: protocol_ID, num_epochs, pre_time, stim_time, tail_time, include_in_analysis
   aligned_response/glom/response  (gloms x frames), attrs mask_vals -> vpn_types.csv names
   acquisition attrs sample_period; stimulus_timing/frame_monitor (channels x 10 kHz samples)
   epochs/epoch_NNN attrs per trial (loom epochs carry 'Loom' in the radius trajectory)

Method (as in their analysis): epoch onsets from the photodiode (frame monitor channel 0
flickers during an epoch, constant between); baseline = mean F in the 1 s pre window;
dF/F = (F - F0)/F0; response amplitude = peak dF/F in the stim window.

Outputs results/turner2022/loom_tuning.json:
   amplitude[glom]  fly-averaged peak dF/F to loom, normalised to the max glom
   gain_sigma       std of log(per-trial population gain) across all loom trials
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
LOOM_SET = ["LC6", "LC26", "LC16", "LPLC2", "LC4", "LPLC1", "LC9", "LC17", "LC12"]   # Turner Fig 5 loom-responsive


def epoch_onsets(fm, fs, n_expected, epoch_s):
    """Onsets (s) of flicker segments in the photodiode trace."""
    win = int(0.05 * fs)
    x = fm[0]
    # moving std via cumulative sums
    c1 = np.cumsum(np.insert(x, 0, 0.0)); c2 = np.cumsum(np.insert(x * x, 0, 0.0))
    m = (c1[win:] - c1[:-win]) / win; v = (c2[win:] - c2[:-win]) / win - m * m
    active = np.sqrt(np.clip(v, 0, None)) > 0.2 * x.std()
    edges = np.flatnonzero(np.diff(active.astype(int)) == 1) + win // 2
    # merge edges closer than half an epoch
    onsets = [edges[0]] if len(edges) else []
    for e in edges[1:]:
        if e - onsets[-1] > 0.5 * epoch_s * fs:
            onsets.append(e)
    return np.array(onsets) / fs


per_fly, per_trial_gain, sources = [], [], []
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
                onsets = epoch_onsets(fm, fs, n_exp, pre + stim + tail)
                if abs(len(onsets) - n_exp) > max(2, 0.05 * n_exp):
                    print(f"  skip {os.path.basename(path)} {fly} {sname} ({proto}): detected {len(onsets)} epochs, expected {n_exp}")
                    continue
                onsets = onsets[:n_exp]
                # which epochs are looms
                ep_keys = sorted(s["epochs"].keys())
                is_loom = []
                for k in ep_keys[:n_exp]:
                    at = " ".join(str(v) for v in s["epochs"][k].attrs.values())
                    is_loom.append(proto == "LoomingSpot" or "Loom" in at)
                t_frames = np.arange(resp.shape[1]) * sp
                peaks = {}
                for on, lo in zip(onsets, is_loom):
                    if not lo:
                        continue
                    base = (t_frames >= on) & (t_frames < on + pre)
                    win = (t_frames >= on + pre) & (t_frames < on + pre + stim)
                    if base.sum() < 3 or win.sum() < 5:
                        continue
                    F0 = resp[:, base].mean(axis=1)
                    dff = (resp[:, win] - F0[:, None]) / np.where(F0 == 0, np.nan, F0)[:, None]
                    pk = np.nanmax(dff, axis=1)
                    for g, v in zip(names, pk):
                        peaks.setdefault(g, []).append(v)
                if not peaks:
                    continue
                df = pd.DataFrame({g: pd.Series(v) for g, v in peaks.items()})
                df = df[[g for g in INCLUDED if g in df.columns]]
                per_fly.append(df.mean(axis=0).rename(f"{os.path.basename(path)}:{fly}:{sname}"))
                loom_cols = [g for g in LOOM_SET if g in df.columns]
                norm = df[loom_cols] / df[loom_cols].mean(axis=0).replace(0, np.nan)
                per_trial_gain.append(norm.mean(axis=1).values)
                sources.append(dict(file=os.path.basename(path), fly=fly, series=sname, protocol=proto,
                                    loom_trials=int(df.shape[0]), gloms=int(df.shape[1])))
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
           series=sources, method="photodiode epoch onsets; dF/F vs 1 s pre; peak in 3 s stim window",
           source="Turner, Krieger, Pang & Clandinin 2022, eLife 11:e82587; Dryad doi:10.5061/dryad.h44j0zpp8",
           extracted_on=datetime.date.today().isoformat())
(OUT / "loom_tuning.json").write_text(json.dumps(out, indent=2))
amp.to_csv(OUT / "loom_peak_dff_by_series.csv")
print("\nnormalised loom amplitude by glomerulus (fly-averaged):")
for k, v in out["amplitude"].items():
    print(f"  {k:<6} {v:.3f}")
print(f"trial gain sigma (log): {out['gain_sigma']}   from {out['gain_n_trials']} loom trials, {out['n_flies']} flies, {out['n_series']} series")
print(f"\nto activate:  cp {OUT/'loom_tuning.json'} data/loom_tuning_measured.json")
