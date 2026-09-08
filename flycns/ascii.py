"""Terminal readouts: circuit sketch, spike raster, regional activity bar."""
import numpy as np
import pandas as pd

CIRCUITS = {
    "gf_escape": """
  loom  ->  LC4/LPLC2 (lobula)  ==chem==>  GF (DNp01)  ==gap 0.8ms==>  TTMn  ->  jump
                                                 \\==gap 0.8ms==>  PSI  ->  DLMn (wing)
  [brain] ........................................ neck ........ [VNC]""",
    "auditory": """
  sound ->  JO-A/JO-B (antenna)  ==gap (calibrated 3 mV)==>  GF (DNp01)   subthreshold alone
  loom  ->  LC4/LPLC2            ==chem=====================>  GF          summates -> spike""",
    "feeding": """
  taste ->  gustatory afferents (MxLbN, aPhN)  ==chem==>  SEZ interneurons (GNG...)  ==chem==>  MN9  ->  proboscis
  wind  ->  JO-FV (control)      must NOT reach MN9""",
}


def sketch(name):
    print(CIRCUITS.get(name, ""))


def raster(model, onsets, readouts, window_ms=120, bin_ms=4, max_trials=3, meta=None):
    """One row per trial per readout type: '.' silent, '|' spike, '#' >1 spike in bin.
    Bottom row: whole-CNS spikes per bin as a 0-9 density digit."""
    sp = model.spike_frame()
    nb = window_ms // bin_ms
    print(f"  raster: {bin_ms} ms/char, {window_ms} ms from onset   (| spike, # multiple)")
    for k, on in enumerate(onsets[:max_trials]):
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + window_ms)].copy()
        w["b"] = ((w.t_ms - on) // bin_ms).astype(int)
        for r in readouts:
            ids = set(model.lif_ids[model.readout_index(r)])
            c = w[w.bodyId.isin(ids)].groupby("b").size()
            row = "".join("#" if c.get(b, 0) > 1 else "|" if c.get(b, 0) == 1 else "." for b in range(nb))
            print(f"  t{k:<2} {r:<6} {row}")
        c = w.groupby("b").size()
        dens = "".join(str(min(9, int(np.log2(c.get(b, 0) + 1)))) if c.get(b, 0) else "." for b in range(nb))
        print(f"  t{k:<2} {'CNS':<6} {dens}   (log2 spikes/bin, 9 = flood)")
        print()


def region_bar(model, onsets, window_ms=250, width=40):
    """Spikes per neuron by superclass during the stimulus windows, as bars."""
    sp = model.spike_frame()
    m = np.zeros(len(sp), bool)
    for on in onsets:
        m |= (sp.t_ms.values >= on) & (sp.t_ms.values < on + window_ms)
    w = sp[m]
    n_by = model.meta.loc[model.lif_ids, "superclass"].value_counts()
    act = (w.groupby("superclass").size() / n_by).dropna().sort_values(ascending=False)
    if not len(act):
        print("  region activity: none"); return
    top = act.max()
    print(f"  region activity (spikes per neuron over {len(onsets)} stimulus windows):")
    for sc, v in act.head(10).items():
        bar = "#" * int(width * v / top) if top > 0 else ""
        print(f"  {sc:<22} {bar:<{width}} {v:.3f}")
