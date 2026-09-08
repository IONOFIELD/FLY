"""
Extract measured loom tuning from Turner, Krieger, Pang & Clandinin 2022 (Dryad
doi:10.5061/dryad.h44j0zpp8) to replace the ordinal LOOM_TUNING placeholder.

Needs (one-time downloads into data/turner2022/):
    datafiles.zip      -> data/turner2022/datafiles/*.hdf5      (22 GB)
    template_brain.zip -> data/turner2022/template_brain/vpn_types.csv
and:  pip install git+https://github.com/ClandininLab/visanalysis git+https://github.com/mhturner/glom_pop

What it computes, following figs/pgs_tuning_pop.py and figs/pgs_trial_variance.py:
    for every ChAT-T2A / Syt1GCaMP6f fly with the PanGlomSuite protocol:
        glom x trial x time dF/F, split by stimulus; the looming-spot stimulus is selected
        by its parameter name ("Loom" in the stimulus/parameter string)
    outputs (results/turner2022/loom_tuning.json):
        amplitude[glom]   mean peak dF/F to loom, fly-averaged, normalised to the max glom
        amplitude_sem[glom]
        gain_sigma        std of log(per-trial population gain), where population gain is
                          the trial's mean normalised amplitude across loom-responsive gloms
        n_flies, n_trials, source, extracted_on
Run:  python fit/extract_turner_loom.py [data/turner2022]
"""
import sys, json, os, glob, datetime
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "data/turner2022")
OUT = Path("results/turner2022"); OUT.mkdir(parents=True, exist_ok=True)
try:
    from visanalysis.analysis import imaging_data, shared_analysis
    import h5py
except ImportError:
    raise SystemExit("pip install git+https://github.com/ClandininLab/visanalysis")

vpn = pd.read_csv(ROOT / "template_brain" / "vpn_types.csv")
val2name = dict(zip(vpn["Unnamed: 0"].values, vpn["vpn_types"].values))
INCLUDED = ["LC4", "LC6", "LC9", "LC11", "LC12", "LC15", "LC16", "LC17", "LC18", "LC21", "LC26", "LPLC1", "LPLC2"]

series = shared_analysis.filterDataFiles(
    data_directory=str(ROOT / "datafiles"),
    target_fly_metadata={"driver_1": "ChAT-T2A", "indicator_1": "Syt1GCaMP6f", "indicator_2": "TdTomato"},
    target_series_metadata={"protocol_ID": "PanGlomSuite", "include_in_analysis": True})
print(f"{len(series)} PanGlomSuite series")

per_fly_amp, per_trial_pop = [], []
loom_label = None
for s in series:
    ID = imaging_data.ImagingDataObject(s["file_name"] + ".hdf5", s["series"], quiet=True)
    with h5py.File(ID.file_path, "r") as f:
        import functools
        from visanalysis.util import h5io
        grp = f.visititems(functools.partial(h5io.find_series, sn=ID.series_number))["aligned_response"]["glom"]
        resp = grp["response"][:]; mask_vals = grp.attrs["mask_vals"]
    names = [val2name.get(int(v), str(v)) for v in mask_vals]
    _, erm = ID.getEpochResponseMatrix(resp, dff=True)                 # (gloms, trials, time)
    uniq, mean_resp, _, trial_by_stim = ID.getTrialAverages(erm)
    # locate the loom stimulus by parameter string
    labels = [str(u) for u in uniq]
    loom_idx = [i for i, l in enumerate(labels) if "loom" in l.lower()]
    if not loom_idx:
        print("  no loom parameter in", labels[:5], "..."); continue
    li = loom_idx[0]; loom_label = labels[li]
    tr = trial_by_stim[li]                                              # (gloms, trials, time)
    peak = np.nanmax(tr, axis=2)                                        # (gloms, trials)
    df = pd.DataFrame(peak.T, columns=names)
    df = df[[n for n in INCLUDED if n in df.columns]]
    per_fly_amp.append(df.mean(axis=0))
    norm = df / df.mean(axis=0).replace(0, np.nan)
    per_trial_pop.append(norm.mean(axis=1).values)                    # per-trial population gain
    print(f"  {os.path.basename(ID.file_path)} s{ID.series_number}: {df.shape[0]} loom trials, gloms {len(df.columns)}")

amp = pd.concat(per_fly_amp, axis=1)                                   # gloms x flies
mean_amp = amp.mean(axis=1); sem_amp = amp.std(axis=1) / np.sqrt(amp.shape[1])
normed = (mean_amp / mean_amp.max()).clip(lower=0)
gains = np.concatenate(per_trial_pop); gains = gains[np.isfinite(gains) & (gains > 0)]
out = dict(amplitude={k: round(float(v), 4) for k, v in normed.items()},
           amplitude_sem={k: round(float(v / mean_amp.max()), 4) for k, v in sem_amp.items()},
           gain_sigma=round(float(np.std(np.log(gains))), 4), gain_n_trials=int(len(gains)),
           n_flies=int(amp.shape[1]), loom_stimulus=loom_label,
           source="Turner, Krieger, Pang & Clandinin 2022, eLife 11:e82587; Dryad doi:10.5061/dryad.h44j0zpp8",
           extracted_on=datetime.date.today().isoformat())
(OUT / "loom_tuning.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
print("\nprotocols.py will now use these measured values (copy to data/loom_tuning_measured.json to activate):")
print(f"  cp {OUT/'loom_tuning.json'} data/loom_tuning_measured.json")
