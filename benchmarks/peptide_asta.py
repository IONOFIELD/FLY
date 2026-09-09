"""
Separating the synaptic and peptidergic effects of Pm3 on the ON pathway (exploratory).

Krieger 2023 (PhD thesis, Stanford, Ch. 2) shows AstA is released by a single visual cell type,
Pm3; AstA-R1 is expressed by L1-L5, C2, Mi1, Mi15, Dm9, Tm2, TmY3 and T2; AstA perfusion or
optogenetic Pm3 activation increases the peak-to-trough dynamic range of Mi1's response to light
flashes; and knocking down AstA-R1 in Mi1 blocks the increase while leaving conventional
transmission intact. In MaleCNS v1.0 Pm3 also inhibits Mi1 with 18,394 GABAergic synapses, so
both channels run between the same cells: the model experiment is the mirror of the genetic one.

DRIVE: Mi1's cholinergic lamina inputs L5 (67,655 synapses) and L3 (34,070), NOT L1. Mi1's
largest input is L1 with 141,873 GLUTAMATERGIC (inhibitory) synapses: the canonical ON response
is a disinhibition (light -> photoreceptors depolarise -> histamine inhibits L1 -> L1 stops
releasing glutamate -> Mi1 released), which a network with no spontaneous activity cannot
express. Driving L3+L5 is the tractable substitute in a silent network and is stated as such.

Arms (mirroring the genetic experiment):
  flash_only               L3+L5 drive; Mi1 baseline response
  flash_plus_Pm3_synaptic  add Pm3 activation, peptide table off (Pm3's GABAergic action alone)
  flash_plus_Pm3_AstA      as above with the declared AstA gain on receptor-expressing types
  AstA_minus_Mi1           as above but Mi1 removed from the target list (AstA-R1 knockdown analogue)

Readout: Mi1 spikes per cell per flash, and the same for the other receptor-expressing types.

OUTCOME (2026-09-09, MaleCNS v1.0): NOT TESTABLE IN THIS MODEL CLASS. Mi1 fires 0.0006 spikes
per cell per flash under 100 Hz drive to its two largest cholinergic inputs, i.e. about one
spike per thousand cells, so it is effectively silent and a multiplicative gain on zero input
is zero. The three modulated arms are therefore identical. The cause is the same one that
blocks the feeding circuit: Mi1 sits under 363,509 inhibitory against 132,350 excitatory
synapses and its real drive is disinhibition from L1, which a network with no spontaneous
activity cannot produce. Pm3 is itself GABAergic, so activating it in a silent network does
nothing even in principle. The declared table and this benchmark are kept as the first thing to
run once cell-specific spontaneous activity exists; the arms are not tuned to make Mi1 fire.

WHAT THIS CANNOT DO: the measured effect is a change in the peak-to-trough shape of a graded,
biphasic calcium response. Every cell in this model spikes and has no graded compartment, so
this is a spike-count analogue of a waveform result, not a reproduction of it. Results are
reported, never scored.

Run:  python benchmarks/peptide_asta.py [n_trials] [data_dir]
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams
from flycns.bench import pulse_protocol, provenance
from flycns.graph import PEPTIDE_MODULATION

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
OUT = Path("results/peptide_asta"); OUT.mkdir(parents=True, exist_ok=True)
DRIVE, HZ = ["L3", "L5"], 100.0     # cholinergic Mi1 inputs; see the note above about L1
PM3_HZ = 50.0                        # optogenetic Pm3 activation analogue
entry = next(e for e in PEPTIDE_MODULATION if e["peptide"] == "AstA")
READOUTS = entry["targets"]

neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")
present = [t for t in READOUTS if (neurons.type == t).any()]
print(f"AstA receptor-expressing types present: {present}")
print(f"Pm3 cells: {int((neurons.type == entry['source']).sum())}; "
      f"Pm3->Mi1 synapses: {int(edges[edges.pre.isin(meta.index[meta.type=='Pm3']) & edges.post.isin(meta.index[meta.type=='Mi1'])].weight.sum()):,}")

def run(label, gain_scale, exclude=(), drive_pm3=False):
    t0 = time.time()
    p = LIFParams(peptide_gain_scale=gain_scale, peptide_exclude=exclude)
    stim = DRIVE + ([entry["source"]] if drive_pm3 else [])
    m = CNSModel(neurons, edges, stim, p, electrical=True)
    rates = {t: HZ for t in DRIVE}
    if drive_pm3:
        rates[entry["source"]] = PM3_HZ
    on = pulse_protocol(m, rates, n_trials=N, pulse_ms=200)
    sp = m.spike_frame()
    out = {}
    for t in present:
        ids = set(m.lif_ids[m.readout_index(t)])
        if not ids:
            continue
        n_cells = len(ids)
        cnt = sum(((sp.t_ms >= o) & (sp.t_ms < o + 300) & sp.bodyId.isin(ids)).sum() for o in on)
        out[t] = round(cnt / n_cells / len(on), 4)
    rec = dict(arm=label, peptide_gain_scale=gain_scale, excluded=list(exclude),
               spikes_per_cell_per_flash=out, pop_rate_hz=float(m.population_rate_hz()),
               peptide_applied=getattr(m, "peptides_applied", []), wall_s=round(time.time() - t0))
    print(f"[{label}] Mi1 {out.get('Mi1', 0):.3f} spikes/cell/flash | "
          f"L1 {out.get('L1', 0):.3f} L2 {out.get('L2', 0):.3f} Tm2 {out.get('Tm2', 0):.3f} | "
          f"CNS {rec['pop_rate_hz']:.4f} Hz | {rec['wall_s']}s")
    return rec

report = dict(benchmark="peptide_asta", drive=DRIVE, drive_hz=HZ, n_trials=N,
              declared=entry, arms={})
report["arms"]["flash_only"] = run("flash_only", 0.0)
report["arms"]["flash_plus_Pm3_synaptic"] = run("flash_plus_Pm3_synaptic", 0.0, drive_pm3=True)
report["arms"]["flash_plus_Pm3_AstA"] = run("flash_plus_Pm3_AstA", 1.0, drive_pm3=True)
report["arms"]["AstA_minus_Mi1"] = run("AstA_minus_Mi1", 1.0, exclude=("Mi1",), drive_pm3=True)

f, a, b, c = (report["arms"][k]["spikes_per_cell_per_flash"] for k in
              ("flash_only", "flash_plus_Pm3_synaptic", "flash_plus_Pm3_AstA", "AstA_minus_Mi1"))
rel = lambda x, base: None if not base else round(x / base, 3)
report["comparison"] = dict(
    Mi1_flash_only=f.get("Mi1"), Mi1_flash_plus_Pm3_synaptic=a.get("Mi1"),
    Mi1_flash_plus_Pm3_AstA=b.get("Mi1"), Mi1_receptor_knockdown=c.get("Mi1"),
    Pm3_synaptic_effect=rel(a.get("Mi1", 0), f.get("Mi1")),
    AstA_effect_over_synaptic=rel(b.get("Mi1", 0), a.get("Mi1")),
    knockdown_returns_to_synaptic=rel(c.get("Mi1", 0), a.get("Mi1")),
    note=("published effect is an increase in the peak-to-trough dynamic range of a graded "
          "biphasic calcium response, blocked by AstA-R1 knockdown in Mi1; this is a spike-count "
          "analogue in a model with no graded compartment, reported and never scored"))
report["provenance"] = provenance()
(OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))
print("\ncomparison:", json.dumps(report["comparison"], indent=1, default=float))
