"""
Graph loading, neurotransmitter signing, literature overrides, null models.

Every deviation from the raw connectome is declared here with a citation.
"""
from pathlib import Path

import numpy as np
import pandas as pd

# Fly neurotransmitter sign convention.
# ACh excitatory; GABA and Glu inhibitory (Liu & Wilson 2013); histamine
# inhibitory (Hardie 1989); amines modulatory and not modelled (sign 0).
SIGN_MAP = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "histamine": -1,
            "dopamine": 0, "octopamine": 0, "serotonin": 0}

# Electrical synapses invisible to EM.
# (pre type, post type, ipsilateral, spikelet_mV, citation); trailing * = type prefix
# spikelet_mV is the postsynaptic jump per presynaptic spike. 9 mV exceeds the
# 7 mV rest-to-threshold gap: a one-to-one relay, as measured for GF->TTMn.
# JON->GF is electrical but subthreshold; hundreds of JONs summate to a
# compound potential of a few mV (Pezier & Blagburn 2013), hence 0.3 mV each.
ELECTRICAL_SYNAPSES = [
    ("DNp01", "TTMn", True, 9.0, "shakB gap junction, 1:1 relay; Tanouye & Wyman 1980; Allen et al. 2006"),
    ("DNp01", "PSI",  True, 9.0, "shakB gap junction, 1:1 relay; Allen et al. 2006; Phelan et al. 2008"),
    # calibrated: population drive at 150 Hz yields a compound GF potential of
    # ~3 mV (Pezier & Blagburn 2013 sound-evoked GF response, subthreshold)
    ("JO-A*", "DNp01", True, {"compound_mV": 3.0, "at_hz": 150}, "JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009"),
    ("JO-B*", "DNp01", True, {"compound_mV": 3.0, "at_hz": 150}, "JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009"),
]
# Mixed synapses that EM annotated as chemical: MaleCNS v1.0 has 679 direct
# JO-A/B -> DNp01 "chemical" synapses. When the electrical model is on, these
# direct edges are removed so the same contact is not counted twice.
MIXED_IN_EM = [("JO-A*", "DNp01"), ("JO-B*", "DNp01")]
TAU_M_S = 0.020   # for compound-potential calibration; must match LIFParams.tau_m_ms
G_GAP_RELAY = 0.2  # ohmic coupling (fraction of leak) for a 1:1 relay pair; population
                   # pairs share this budget so many resting partners cannot shunt the post cell

# Gustatory afferents by ANATOMY (MaleCNS v1.0 `subclass`); modality is unannotated.
GUSTATORY_SUBCLASSES = ["labellar bristle", "taste peg", "pharyngeal sensillum"]

# Regional gain: synapses from SEZ-intrinsic neurons get their own multiplier.
# Definition (since 2026-09-08): more than half of a neuron's synapses lie in SEZ
# compartments (GNG, PRW, SAD, FLA, CAN, AMMC, PENP), from the connectome's own
# roiInfo (fit/roi_membership.py -> data/roi_membership.parquet), AND the neuron is
# not sensory, motor or efferent (those are inputs to and outputs of the region, not
# local processing). 3,519 neurons, 4.0% of synaptic weight. The earlier definition
# was a type-name prefix rule (2,280 neurons, 3.6%); it is kept as a fallback when
# roi_membership.parquet is absent, and the two overlap in 1,414 neurons.
# Reading: SEZ local synapses are effectively stronger than synapse counts
# imply (hypothesis to test in vivo), needed for taste -> MN9 propagation.
SEZ_TYPE_PREFIXES = ("GNG", "SAD", "PRW", "FLA", "CAN")     # fallback only
SEZ_ROI_PREFIXES = ("GNG", "PRW", "SAD", "FLA", "CAN", "AMMC", "PENP")
_SEZ_BODIES = None


def sez_bodies(neurons=None, data_dir="data"):
    """Body IDs of SEZ-intrinsic neurons by anatomy; empty if the membership file is absent."""
    global _SEZ_BODIES
    if _SEZ_BODIES is None:
        f = Path(data_dir) / "roi_membership.parquet"
        if not f.exists():
            _SEZ_BODIES = frozenset()
        else:
            m = pd.read_parquet(f).set_index("bodyId")
            keep = set(m.index[m.sez_frac > 0.5])
            if neurons is not None:
                sc = neurons.set_index("bodyId")["superclass"].fillna("")
                keep &= set(sc.index[~sc.str.contains("sensory|motor|efferent")])
            _SEZ_BODIES = frozenset(keep)
    return _SEZ_BODIES
# Regional gain is DISABLED (1.0) since 2026-09-08. It was fitted to 2.0 under a
# type-name prefix definition of the SEZ, where it appeared to open the taste -> MN9
# pathway. With the SEZ defined anatomically from the connectome's own compartment
# annotations, the population is net inhibitory onto that pathway under either
# definition, so the gain amplifies inhibition along with excitation and MN9 stays
# silent at every value (benchmarks/fit_regional.py). The underlying reason is that
# the pathway is disinhibitory (benchmarks/feeding_mechanism.py); no uniform
# manipulation reproduces it. Kept as a parameter for anyone who wants to sweep it.
SEZ_GAIN_FITTED = 1.0

# Screen-selected gustatory subset (benchmarks/screen_afferents.py at w_scale 1.0,
# Sept 2026): the gustatory types that individually drive MN9. This is the closest
# proxy for sugar GRNs the annotation permits; it is selected by wiring, NOT a
# modality annotation, and is labelled as such wherever it is used.
GUSTATORY_MN9_DRIVING = ["aPhM2a", "aPhM5", "PhG1c", "claw_tpGRN", "LB3a", "LB3b", "LB3c", "LB3d"]

# Known feeding circuit types (Shiu et al. 2022 eLife; Shiu et al. 2024 Nature).
# MaleCNS naming may differ; benchmarks/feeding.py discovers what is present.
FEEDING_SECOND_ORDER = ["Fdg", "Clavicle", "Zorro", "Rattle", "Phantom", "Bract",
                        "Roundup", "Usnea", "Cleaver"]
FEEDING_MOTOR = ["MN9", "MN6", "MN8", "MN11"]
SUGAR_GRN_PATTERNS = [r"Gr64f", r"sugar", r"GRN.*sugar", r"sugar.*GRN", r"Gr5a", r"Taste.*sugar", r"sweet"]
BITTER_GRN_PATTERNS = [r"Gr66a", r"bitter", r"GRN.*bitter", r"Taste.*bitter"]
TASTE_SENSORY_FALLBACK = ["BM_Taste"]   # MaleCNS labellar taste population, modality unsplit

# Per-type intrinsic overrides. (type, parameter, value, citation)
# GF fires 1 to 2 spikes per loom regardless of input strength
# (von Reyn et al. 2014; Ache et al. 2019): spike-triggered adaptation.
INTRINSIC_OVERRIDES = [
    ("DNp01", "b_adapt_mV", 15.0, "all-or-none GF response; von Reyn 2014; Ache 2019. Fitted by pre-stated rule "
                                  "(smallest value giving 1-2 spikes/response, <=2 on strong looms, P(response) 0.5-1, "
                                  "TTMn relay intact; benchmarks/fit_gf_adapt.py). 30 mV under the ordinal loom tuning "
                                  "(2026-09-08), 15 mV under the measured Turner 2022 tuning (2026-09-08); 0 mV gives up to 8 spikes."),
    # superclass-level: sustained fly motor neuron firing stays well under 100 Hz
    # (leg MNs, Azevedo et al. 2020; proboscis MNs, McKellar et al. 2020). Without
    # adaptation the LIF motor pool pins at the 450 Hz refractory ceiling.
    ("superclass:cb_motor", "b_adapt_mV", 10.0, "sustained MN rates < 100 Hz; Azevedo 2020; McKellar 2020"),
    # vnc_motor deliberately NOT adapted: TTMn/DLMn follow GF spikes 1:1 up to 100 Hz
    # (Tanouye & Wyman 1980); an adaptation term there blocked relay of GF doublets
    # once measured loom tuning broadened the drive (2026-09-08).
]


def load_graph(data_dir="data", min_weight=5):
    import os
    neurons = pd.read_parquet(f"{data_dir}/neurons.parquet")
    edges = pd.read_parquet(f"{data_dir}/edges.parquet")
    edges = edges[edges["weight"] >= min_weight].copy()
    edges["sign"] = edges["nt"].str.lower().map(SIGN_MAP)
    n_unknown = edges["sign"].isna().sum()
    edges = edges[edges["sign"].notna() & (edges["sign"] != 0)].reset_index(drop=True)
    edges["sign"] = edges["sign"].astype(float)
    # robustness: flip a random fraction of neurotransmitter signs (FLYCNS_SIGNFLIP, e.g. 0.05)
    flip = float(os.environ.get("FLYCNS_SIGNFLIP", "0"))
    if flip > 0:
        rng = np.random.default_rng(int(os.environ.get("FLYCNS_SEED", "0")) + 1000)
        m = rng.random(len(edges)) < flip
        edges.loc[m, "sign"] *= -1
        print(f"[graph] SIGN-FLIP PERTURBATION: {m.sum():,} edges ({flip:.0%}) sign-inverted")
    print(f"[graph] {len(neurons):,} neurons, {len(edges):,} signed edges "
          f"(dropped {n_unknown:,} unknown-NT pairs and modulatory amines)")
    return neurons, edges


def pre_class(superclass):
    """Three presynaptic classes for class-wise synaptic gain (regional fitting)."""
    sc = "" if superclass is None or superclass != superclass else str(superclass)
    if "sensory" in sc:
        return "sensory"
    if "descending" in sc or "ascending" in sc:
        return "relay"
    return "local"


def type_region(t):
    """Fallback region tag from the type-name prefix (used only without roi_membership.parquet)."""
    t = "" if t is None or t != t else str(t)
    return "SEZ" if t.startswith(SEZ_TYPE_PREFIXES) else "other"


def region_of(neurons, data_dir="data"):
    """Per-bodyId region ('SEZ'/'other'), anatomical where available, else prefix rule."""
    bodies = sez_bodies(neurons, data_dir)
    idx = neurons.set_index("bodyId")
    if bodies:
        r = pd.Series("other", index=idx.index)
        r.loc[list(bodies & set(idx.index))] = "SEZ"
        return r, f"anatomical (roiInfo, {len(bodies):,} SEZ-intrinsic neurons)"
    return idx["type"].map(type_region), "type-name prefix fallback"


def gustatory_afferents(neurons):
    """Types whose subclass is anatomically gustatory (requires fetch_annotations.py)."""
    if "subclass" not in neurons:
        return []
    m = neurons[neurons["subclass"].isin(GUSTATORY_SUBCLASSES)]
    return sorted(m["type"].dropna().unique())


def infer_side(meta):
    """somaSide, filled from instance suffix (_L/_R) then soma x relative to the
    midline of labelled neurons. Sensory neurons often lack soma coordinates."""
    side = meta["somaSide"].copy()
    inst = meta["instance"].fillna("")
    side = side.fillna(inst.str.extract(r"_([LR])(?:$|_|\b)")[0])
    if "x" in meta and side.isna().any():
        lab = meta[side.isin(["L", "R"]) & meta["x"].notna()]
        if len(lab):
            mid = lab.groupby(side[lab.index])["x"].median().mean()
            left_is_low = lab.loc[side[lab.index] == "L", "x"].median() < mid
            guess = np.where(meta["x"] < mid, "L" if left_is_low else "R", "R" if left_is_low else "L")
            side = side.fillna(pd.Series(guess, index=meta.index).where(meta["x"].notna()))
    return side


def pre_class(superclass):
    """Three presynaptic classes for class-wise synaptic gain (regional fitting)."""
    sc = "" if superclass is None or superclass != superclass else str(superclass)
    if "sensory" in sc:
        return "sensory"
    if "descending" in sc or "ascending" in sc:
        return "relay"
    return "local"


def type_region(t):
    """Fallback region tag from the type-name prefix (used only without roi_membership.parquet)."""
    t = "" if t is None or t != t else str(t)
    return "SEZ" if t.startswith(SEZ_TYPE_PREFIXES) else "other"


def region_of(neurons, data_dir="data"):
    """Per-bodyId region ('SEZ'/'other'), anatomical where available, else prefix rule."""
    bodies = sez_bodies(neurons, data_dir)
    idx = neurons.set_index("bodyId")
    if bodies:
        r = pd.Series("other", index=idx.index)
        r.loc[list(bodies & set(idx.index))] = "SEZ"
        return r, f"anatomical (roiInfo, {len(bodies):,} SEZ-intrinsic neurons)"
    return idx["type"].map(type_region), "type-name prefix fallback"


def gustatory_afferents(neurons):
    """Types whose subclass is anatomically gustatory (requires fetch_annotations.py)."""
    if "subclass" not in neurons:
        return []
    m = neurons[neurons["subclass"].isin(GUSTATORY_SUBCLASSES)]
    return sorted(m["type"].dropna().unique())


def infer_side(meta):
    """somaSide where annotated; else the instance suffix (_L/_R, which MaleCNS
    gives even to sensory neurons); else soma x against the midline; else None."""
    side = meta["somaSide"].copy()
    if "instance" in meta:
        suf = meta["instance"].fillna("").str.extract(r"_([LR])$")[0]
        side = side.where(side.isin(["L", "R"]), suf)
    if "x" in meta:
        mid = meta.loc[meta["somaSide"].isin(["L", "R"]), "x"].median()
        guess = pd.Series(np.where(meta["x"] < mid, "L", "R"), index=meta.index)
        guess[meta["x"].isna()] = None
        side = side.where(side.isin(["L", "R"]), guess)
    return side


def _match(meta, pattern):
    t = meta["type"].fillna("")
    return meta[t.str.startswith(pattern[:-1])] if pattern.endswith("*") else meta[t == pattern]


def electrical_pairs(neurons):
    meta = neurons.set_index("bodyId").copy()
    meta["somaSide"] = infer_side(meta)
    pairs = []
    # pool prefix entries with the same post so calibration counts the whole population
    pooled = {}
    for pre_t, post_t, ipsi, spk, src in ELECTRICAL_SYNAPSES:
        key = (post_t, ipsi, str(spk))
        pooled.setdefault(key, []).append((pre_t, spk, src))
    for (post_t, ipsi, _), entries in pooled.items():
        pres = pd.concat([_match(meta, e[0]) for e in entries])
        spk, src = entries[0][1], "; ".join(sorted({e[2] for e in entries}))
        posts = _match(meta, post_t)
        g_gap = G_GAP_RELAY
        if isinstance(spk, dict):
            # steady-state summed depolarization = n_pre * rate * spikelet * tau_m.
            # Count the neurons that will actually converge on ONE post cell.
            sided = pres["somaSide"].isin(["L", "R"])
            if ipsi and sided.any():
                n_conv = pres[sided].groupby("somaSide").size().mean() + (~sided).sum()
            else:
                n_conv = len(pres)
            n_conv = max(float(n_conv), 1.0)
            spk = spk["compound_mV"] / (n_conv * spk["at_hz"] * TAU_M_S)
            g_gap = G_GAP_RELAY / n_conv
        assert np.isfinite(spk) and np.isfinite(g_gap), f"bad electrical calibration for {post_t}"
        for b, row in pres.iterrows():
            cand = posts[posts["somaSide"] == row["somaSide"]] if ipsi else posts
            if not len(cand):
                cand = posts
            for q in cand.index:
                pairs.append(dict(pre=b, post=q, pre_type=row["type"], post_type=post_t,
                                  spikelet_mV=float(spk), g_gap=float(g_gap), source=src))
    return pd.DataFrame(pairs)


def drop_mixed_chemical(edges, neurons):
    """Remove direct chemical edges for pairs declared in MIXED_IN_EM."""
    meta = neurons.set_index("bodyId")
    keep = pd.Series(True, index=edges.index)
    n_drop = 0
    for pre_t, post_t in MIXED_IN_EM:
        pres = _match(meta, pre_t).index
        posts = _match(meta, post_t).index
        m = edges["pre"].isin(pres) & edges["post"].isin(posts)
        n_drop += int(edges.loc[m, "weight"].sum())
        keep &= ~m
    if n_drop:
        print(f"[graph] removed {n_drop:,} EM chemical synapses on declared mixed pairs (modelled as electrical)")
    return edges[keep].reset_index(drop=True)


def rewire_null(edges, seed=0):
    """Degree-preserving null: permute postsynaptic targets, keep out-degree,
    synapse counts and NT sign of every presynaptic neuron."""
    rng = np.random.default_rng(seed)
    out = edges.copy()
    out["post"] = rng.permutation(out["post"].values)
    return out[out["pre"] != out["post"]].reset_index(drop=True)
