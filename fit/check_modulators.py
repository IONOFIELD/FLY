"""
Read-only checks against MaleCNS v1.0. Safe to run while a benchmark battery is going:
reads parquet only, writes results/checks/ and nothing else.

PART A. Do the candidate state-signalling neurons from the literature exist here by name?
  IPCs (insulin-producing cells, pars intercerebralis)      eLife 2024;13:e98514
  SEZ octopaminergic AKHR cells                             eLife 2016;5:e15693
  hugin neurons (SEZ)                                       Melcher & Pankratz 2005
  leucokinin, NPF, sNPF, AstA, DSK, AKH
  SLC5A11 / CN nutrient sensor (ellipsoid body)
A "likely yes" from anatomical inference is not good enough; this answers by name.

PART B. Is there a TH-VUM-shaped body among the untyped GNG cells?
Marella, Mann & Scott 2012 (Neuron 73:941-950) describe TH-VUM as: a SINGLE ventral
unpaired median interneuron, soma on the midline in the SEZ, with BROAD arborisation
across the gnathal ganglion, and nSyb-GFP labelling ALL arbors, i.e. almost entirely
PRESYNAPTIC. That is a searchable morphological signature:
  - no somaSide (midline), or an instance suffix _M
  - >50% of synapses in SEZ compartments
  - high pre:post ratio (output-dominated)
  - large total synapse count (broad arbor), few or one such cell
Candidates are ranked, not identified. Confirmation would need NBLAST against the
light-level template, which this script cannot do.

Run:  python fit/check_modulators.py [data_dir]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "data")
OUT = Path("results/checks")
OUT.mkdir(parents=True, exist_ok=True)

neurons = pd.read_parquet(DATA / "neurons.parquet")
meta = neurons.set_index("bodyId")
edges = pd.read_parquet(DATA / "edges.parquet")
roi_f = DATA / "roi_membership.parquet"
roi = pd.read_parquet(roi_f).set_index("bodyId") if roi_f.exists() else None
if roi is None:
    raise SystemExit("run fit/roi_membership.py first (needs data/roi_membership.parquet)")

t = neurons["type"].fillna("")
inst = neurons["instance"].fillna("") if "instance" in neurons else pd.Series("", index=neurons.index)
report = {}

# ---------------------------------------------------------------- PART A
print("=" * 70)
print("PART A: do the literature's state-signalling neurons exist in MaleCNS by name?")
print("=" * 70)

PATTERNS = {
    "IPC / insulin-producing": [r"^IPC", r"insulin", r"DILP", r"Ilp"],
    "AKH / adipokinetic": [r"AKH"],
    "octopaminergic VUM (SEZ-adjacent)": [r"^OA-VUM", r"^OA-VPM", r"^OA-AL"],
    "hugin": [r"hugin", r"^HUG"],
    "leucokinin": [r"^LHLK", r"^SELK", r"leucokinin", r"^LK$"],
    "NPF / sNPF": [r"NPF", r"sNPF"],
    "allatostatin A": [r"AstA", r"allatostatin"],
    "drosulfakinin": [r"DSK", r"drosulfakinin"],
    "SLC5A11 / CN nutrient sensor": [r"SLC5A11", r"^CN$", r"^ER[0-9]", r"^EPG$"],
    "clock (reference, known present)": [r"LNv", r"LNd", r"^DN1", r"^DN2", r"^DN3"],
}
rows = []
for label, pats in PATTERNS.items():
    hit = pd.Series(False, index=neurons.index)
    for p in pats:
        hit |= t.str.contains(p, case=False, regex=True) | inst.str.contains(p, case=False, regex=True)
    types = sorted(neurons.loc[hit, "type"].dropna().unique())
    n_cells = int(hit.sum())
    rows.append(dict(group=label, n_cells=n_cells, n_types=len(types), types=types[:8]))
    status = "FOUND" if n_cells else "ABSENT"
    print(f"  {status:<7} {label:<34} {n_cells:>5} cells, {len(types):>3} types  {types[:5]}")
report["part_a_by_name"] = rows

# transmitter route: monoamine cells that are actually present, by region
print("\n  monoamine cells by transmitter and whether any sit in the SEZ:")
mono = neurons[neurons["consensusNt"].isin(["dopamine", "octopamine", "serotonin"])].set_index("bodyId")
mono = mono.join(roi[["sez_frac", "dominant_roi"]])
mono_tab = []
for nt in ["dopamine", "octopamine", "serotonin"]:
    sub = mono[mono["consensusNt"] == nt]
    in_sez = int((sub["sez_frac"] > 0.5).sum())
    near = int((sub["sez_frac"] > 0.2).sum())
    mono_tab.append(dict(transmitter=nt, n=len(sub), sez_frac_gt_0_5=in_sez, sez_frac_gt_0_2=near,
                         top_rois=sub["dominant_roi"].value_counts().head(3).to_dict()))
    print(f"    {nt:<12} {len(sub):>5} cells, >50% SEZ: {in_sez:>3}, >20% SEZ: {near:>3}, "
          f"top ROIs {list(sub['dominant_roi'].value_counts().head(3).index)}")
report["part_a_monoamines"] = mono_tab

# ---------------------------------------------------------------- PART B
print("\n" + "=" * 70)
print("PART B: TH-VUM-shaped candidates among untyped GNG bodies")
print("=" * 70)

j = meta.join(roi[["sez_frac", "dominant_roi", "n_syn"]])
pre = edges.groupby("pre")["weight"].sum().rename("out_w")
post = edges.groupby("post")["weight"].sum().rename("in_w")
j = j.join(pre).join(post)
j["out_w"] = j["out_w"].fillna(0)
j["in_w"] = j["in_w"].fillna(0)
j["total_w"] = j["out_w"] + j["in_w"]
j["pre_frac"] = j["out_w"] / j["total_w"].replace(0, np.nan)

side = j["somaSide"]
suffix_m = j["instance"].fillna("").str.endswith("_M")
midline = side.isna() | suffix_m

cand = j[(j["sez_frac"] > 0.5) & midline & (j["total_w"] >= 200)].copy()
cand = cand.sort_values(["pre_frac", "total_w"], ascending=False)
print(f"  midline bodies with >50% SEZ synapses and >=200 synapses: {len(cand)}")

strict = cand[(cand["pre_frac"] > 0.7)]
print(f"  of those, output-dominated (pre_frac > 0.7): {len(strict)}")
cols = ["type", "instance", "consensusNt", "predictedNt", "dominant_roi", "sez_frac",
        "out_w", "in_w", "pre_frac"]
show = strict[cols].head(20).copy()
show["sez_frac"] = show["sez_frac"].round(2)
show["pre_frac"] = show["pre_frac"].round(2)
print(show.to_string())

# which of these are dopamine-predicted at all
dop = strict[(strict["predictedNt"] == "dopamine") | (strict["consensusNt"] == "dopamine")]
print(f"\n  of the output-dominated midline SEZ bodies, dopamine-predicted: {len(dop)}")
if len(dop):
    print(dop[cols].head(10).to_string())

# for the top candidates, do they contact MN9's inhibitory relays? (the disinhibitory arm)
mn9 = meta.index[meta["type"] == "MN9"]
sign_map = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "histamine": -1}
e = edges.assign(sgn=edges["nt"].str.lower().map(sign_map).fillna(0))
inh_relays = set(e[(e.post.isin(mn9)) & (e.sgn < 0)].pre)
exc_relays = set(e[(e.post.isin(mn9)) & (e.sgn > 0)].pre)
top = list(strict.index[:40])
onto = e[e.pre.isin(top)]
hits = []
for b in top:
    d = onto[onto.pre == b]
    w_inh = int(d[d.post.isin(inh_relays)].weight.sum())
    w_exc = int(d[d.post.isin(exc_relays)].weight.sum())
    if w_inh or w_exc:
        hits.append(dict(bodyId=int(b), type=strict.at[b, "type"], pre_frac=round(float(strict.at[b, "pre_frac"]), 2),
                         onto_MN9_inhibitory=w_inh, onto_MN9_excitatory=w_exc))
hits = sorted(hits, key=lambda r: -(r["onto_MN9_inhibitory"] + r["onto_MN9_excitatory"]))
print(f"\n  candidates that contact MN9's relays at all: {len(hits)}")
for h in hits[:10]:
    print(f"    body {h['bodyId']} type {str(h['type']):<10} pre_frac {h['pre_frac']}  "
          f"-> MN9-inhibitory {h['onto_MN9_inhibitory']:>6}  -> MN9-excitatory {h['onto_MN9_excitatory']:>6}")

report["part_b"] = dict(
    n_midline_sez=int(len(cand)), n_output_dominated=int(len(strict)),
    n_dopamine_predicted=int(len(dop)),
    top_candidates=[{k: (float(v) if isinstance(v, (np.floating, float)) else
                         int(v) if isinstance(v, (np.integer,)) else v)
                     for k, v in r.items()}
                    for r in strict[cols].head(20).reset_index().to_dict("records")],
    contacts_mn9_relays=hits[:20],
    signature=("Marella 2012: single midline VUM, broad GNG arbor, nSyb-GFP labels all arbors "
               "(output-dominated). Candidates are ranked by that signature only; confirming "
               "identity needs NBLAST against the light-level template."))

(OUT / "modulator_checks.json").write_text(json.dumps(report, indent=2, default=str))
print(f"\nwrote {OUT / 'modulator_checks.json'}")
