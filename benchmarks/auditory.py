"""
Benchmark 3: Johnston's organ -> giant fiber (mechanosensory input to escape).

Why this circuit: JON -> GF is a documented MIXED synapse that is PRIMARILY
ELECTRICAL (shakB); shakB2 mutants lose the JON-evoked GF current (Pezier &
Blagburn 2013; Yorozu et al. 2009). In vivo, sound alone leaves GF
subthreshold and summates with visual drive (von Reyn 2014). That yields
ablation predictions with opposite signs to the GF->TTMn case:

  A1 JO-A/B drive alone (pulse-song rates): GF response probability low, <= 0.3
  A2 JO drive + loom: GF latency shorter, or hit rate/count higher, than loom alone
  A3 pathway ablation: MaleCNS EM annotates the mixed JON-GF contact as
     chemical (679 synapses). Removing BOTH the electrical model and those
     EM edges must abolish the JO-evoked GF depolarisation (measured from GF
     membrane voltage, Pezier & Blagburn 2013); with the electrical model it
     must exceed 0.5 mV. "chem_only_EM" (EM as annotated, no electrical) is
     reported for reference: it double-counts nothing but mislabels the contact.
  A4 rewired null, chemical only: no JO -> GF depolarisation (the electrical
     pairs are declared anatomy and bypass rewiring, so they are removed here)
  A5 CNS stays quiet (< 0.05 Hz/neuron)

JO-A and JO-B are the vibration-sensitive subgroups (Kamikouchi et al. 2009);
pulse song drives them at ~100-200 Hz burst rates.
Writes results/auditory/{report.json, provenance.md, *.csv}
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, rewire_null, CNSModel, LIFParams, loom_protocol
from flycns.graph import drop_mixed_chemical
from flycns.protocols import LOOM_TUNING
from flycns.bench import find_types, pulse_protocol
from flycns import ascii as A

OUT = Path("results/auditory"); OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
JO_HZ = 150

neurons, edges = load_graph(DATA)
A.sketch("auditory")
meta = neurons.set_index("bodyId")
jo_types = [t for t in find_types(neurons, [r"^JO-A", r"^JO-B"]) if not t.endswith("unclear")]
print("JO vibration types:", jo_types)
if not jo_types:
    raise SystemExit("no JO-A/JO-B types found")
report = {"benchmark": "auditory", "n_trials": N_TRIALS, "jo_types": jo_types, "jo_hz": JO_HZ,
          "arms": {}, "checks": {}}


def gf_stats(model, onsets, window_ms=160):
    sp = model.spike_frame()
    gf_ids = set(model.lif_ids[model.readout_index("DNp01")])
    rows = []
    for on in onsets:
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + window_ms) & sp.bodyId.isin(gf_ids)]
        rows.append(dict(gf=len(w) / max(len(gf_ids), 1), hit=int(len(w) > 0),
                         lat=(w.t_ms.min() - on) if len(w) else np.nan))
    df = pd.DataFrame(rows)
    return dict(gf_per_cell=float(df.gf.mean()), gf_hit=float(df.hit.mean()),
                gf_lat_ms=float(df.lat.mean()), pop_rate_hz=float(model.population_rate_hz()))


def run(label, stim_types, electrical=True, edge_df=edges, jo=True, loom=False):
    t0 = time.time()
    m = CNSModel(neurons, edge_df, stim_types, LIFParams(), electrical=electrical,
                 monitor_types=("DNp01",))
    if loom and jo:
        # loom protocol drives LC types; add JO at fixed rate during each loom burst
        tuning = m.stim["type"].map(LOOM_TUNING).fillna(0).values
        is_jo = m.stim["type"].isin(jo_types).values
        rng = np.random.default_rng(0); onsets = []
        for _ in range(N_TRIALS):
            m.run(300); onsets.append(m.t_ms)
            gain = rng.lognormal(0, 0.5)
            m.set_stim_rates(tuning * 5.0 * gain + is_jo * JO_HZ); m.run(60)
            m.set_stim_rates(np.zeros(len(m.stim)))
        m.run(300)
    elif loom:
        onsets = [on for on, _ in loom_protocol(m, n_trials=N_TRIALS)]
    else:
        onsets = pulse_protocol(m, {t: JO_HZ for t in stim_types}, n_trials=N_TRIALS, pulse_ms=60)
    s = gf_stats(m, onsets); s["wall_s"] = round(time.time() - t0, 1)
    A.raster(m, onsets, ["DNp01"], window_ms=100, max_trials=2)
    s["gf_depol_mV"] = float(np.nanmean(m.peak_depolarization_mV(onsets, 100)))
    print(f"[{label}] GF {s['gf_per_cell']:.2f}/cell hit {s['gf_hit']:.2f} lat {s['gf_lat_ms']:.1f} ms "
          f"depol {s['gf_depol_mV']:.2f} mV | CNS {s['pop_rate_hz']:.4f} Hz | {s['wall_s']}s")
    report["arms"][label] = s
    return s


a1 = run("A1_jo_alone_electrical", jo_types, electrical=True)
a3 = run("A3_jo_alone_no_electrical_no_mixed_edges", jo_types, electrical=False,
         edge_df=drop_mixed_chemical(edges, neurons))
ref = run("ref_jo_alone_chem_only_EM_as_annotated", jo_types, electrical=False)
lo = run("loom_alone", list(LOOM_TUNING), electrical=True, jo=False, loom=True)
a2 = run("A2_loom_plus_jo", list(LOOM_TUNING) + jo_types, electrical=True, jo=True, loom=True)
a4 = run("A4_jo_rewired_null_chem", jo_types, electrical=False,
         edge_df=rewire_null(drop_mixed_chemical(edges, neurons)))

report["checks"] = {
    "A1_jo_alone_mostly_subthreshold": a1["gf_hit"] <= 0.3,
    "A2_summation_with_loom": (a2["gf_hit"] > lo["gf_hit"]) or (a2["gf_per_cell"] > lo["gf_per_cell"])
                              or (a2["gf_lat_ms"] < lo["gf_lat_ms"] - 1.0),
    "A3_declared_pairs_carry_jo_drive": (a1["gf_depol_mV"] > 0.5) and (a3["gf_depol_mV"] < 0.5 * a1["gf_depol_mV"]),
    "A4_null_silent": a4["gf_hit"] < 0.1 and a4["gf_depol_mV"] < 0.5,
    "A5_cns_quiet": max(a1["pop_rate_hz"], a2["pop_rate_hz"]) < 0.05,
}
report["pass"] = all(report["checks"].values())
(OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))
prov = ["# Provenance: auditory (JON -> GF)", "", f"JO types: {jo_types} at {JO_HZ} Hz",
        "Electrical JON->GF calibrated to a 3 mV compound GF potential at 150 Hz (Pezier & Blagburn 2013)",
        "", "## Checks"] + [f"- {k}: {'PASS' if v else 'FAIL'}" for k, v in report["checks"].items()]
(OUT / "provenance.md").write_text("\n".join(prov))
print("\nCHECKS"); [print(f"  {'PASS' if v else 'FAIL'}  {k}") for k, v in report["checks"].items()]
print(f"\nBENCHMARK {'PASS' if report['pass'] else 'FAIL'}  -> {OUT}/report.json")
