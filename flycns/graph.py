"""
Graph loading, neurotransmitter signing, literature overrides, null models.

Every deviation from the raw connectome is declared here with a citation.
"""
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
    ("DNp01", "b_adapt_mV", 30.0, "all-or-none GF response; von Reyn 2014; Ache 2019"),
]


def load_graph(data_dir="data", min_weight=5):
    neurons = pd.read_parquet(f"{data_dir}/neurons.parquet")
    edges = pd.read_parquet(f"{data_dir}/edges.parquet")
    edges = edges[edges["weight"] >= min_weight].copy()
    edges["sign"] = edges["nt"].str.lower().map(SIGN_MAP)
    n_unknown = edges["sign"].isna().sum()
    edges = edges[edges["sign"].notna() & (edges["sign"] != 0)].reset_index(drop=True)
    edges["sign"] = edges["sign"].astype(float)
    print(f"[graph] {len(neurons):,} neurons, {len(edges):,} signed edges "
          f"(dropped {n_unknown:,} unknown-NT pairs and modulatory amines)")
    return neurons, edges


def pre_class(superclass):
    """Three presynaptic classes for class-wise synaptic gain (tier 2 fitting)."""
    sc = "" if superclass is None or superclass != superclass else str(superclass)
    if "sensory" in sc:
        return "sensory"
    if "descending" in sc or "ascending" in sc:
        return "relay"
    return "local"


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
    """Three presynaptic classes for class-wise synaptic gain (tier 2 fitting)."""
    sc = "" if superclass is None or superclass != superclass else str(superclass)
    if "sensory" in sc:
        return "sensory"
    if "descending" in sc or "ascending" in sc:
        return "relay"
    return "local"


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
