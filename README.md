# flycns: connectome-constrained LIF benchmarks on the Drosophila MaleCNS

One-command, repeatable simulation of the Janelia FlyEM male CNS connectome
(v1.0, 166,700 neurons, brain + ventral nerve cord) as a signed leaky
integrate-and-fire network, validated against published circuit physiology.

## Quick start
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    export NEUPRINT_TOKEN="..."         # neuprint.janelia.org account page
    python fetch_malecns.py              # ~5 min, writes data/*.parquet
    ./run_all.sh 20                      # benchmark + checks + 3D cascade

## Layout
    flycns/graph.py      load, NT signing, electrical synapses, intrinsic overrides, null model
    flycns/model.py      Brian2 LIF network (chemical + electrical synapses, adaptation)
    flycns/protocols.py  stimulus protocols with citations (loom after Turner, Krieger et al. 2022)
    benchmarks/          one file per circuit; writes results/<name>/report.json + provenance.md
    viz/cascade_3d.py    animated 3D cascade from saved spikes (plotly, html)
    make_testdata.py     synthetic graph for smoke tests (no data required)

## Benchmarks
| circuit | stimulus | readout | criteria | status |
|---|---|---|---|---|
| gf_escape | loom, 9 LC types (Turner et al. 2022) | DNp01, TTMn, PSI | 1-2 GF spikes per response, hit >= 0.8, 1:1 TTMn relay at 0.8 ms, chemical-only ablation kills relay, rewired null silent, quiet CNS | 8 checks |
| auditory | JO-A/B at 150 Hz | DNp01 membrane V | JO alone subthreshold (~3 mV, Pezier & Blagburn 2013), summates with loom, electrical carries it, chem-only+rewired null silent | 5 checks |
| feeding | sugar GRN (or Fdg proxy), bitter | MN9 by side | contra > ipsi, monotonic 10-200 Hz, bitter suppresses, second-order sufficient, null silent (port of Shiu 2024) | 5 checks |

Feeding status on v1.0: under a single global synaptic scale the model does not reproduce
sugar-driven MN9 activation. `benchmarks/screen_afferents.py` (a sufficiency screen after
Shiu 2024) shows no afferent type drives MN9 at w_scale 0.3, and at Shiu's 1.0 MN9 responds
to auditory and vibration afferents as strongly as to taste while the CNS floods. Taste
afferents suppress MN9 when co-applied at intermediate gain. `benchmarks/fit_gains.py`
searches class-wise gains (sensory, relay, local) against all circuits at once.

`benchmarks/summarize.py` writes `results/SUITE.md`: every check across every circuit under one declared parameter set.

## Every deviation from the raw connectome is declared
See `flycns/graph.py` (SIGN_MAP, ELECTRICAL_SYNAPSES, INTRINSIC_OVERRIDES) and
`results/*/provenance.md`. Effective chemical kick is 0.275 mV (Shiu 2024 unitary
weight) x 0.3 (MaleCNS rescale from sweep) = 0.0825 mV per synapse; chemical delay
1.8 ms (Shiu 2024). Electrical synapses: 1:1 relays use a 9 mV spikelet with 0.8 ms
latency (GF-TTMn); population inputs (JON-GF) are calibrated to a measured compound
potential and share one coupling budget so resting partners cannot shunt the target.
Monoamines are sign 0 by deliberate scope (Shiu used +1; no principled default exists).

## What this builds on and what it adds
This is a port and extension of Shiu et al. 2024 (Nature) from FlyWire to MaleCNS.
The model, its parameters, and the feeding validation design are theirs; the giant
fiber physiology is from von Reyn, Ache, Tanouye & Wyman, Allen, and Augustin; the
loom stimulus follows Turner, Krieger, Pang & Clandinin 2022. What this repository
adds is a test harness: explicit pass/fail criteria per circuit, ablation and null
arms, one shared parameter set, and a declared list of every deviation from the raw
connectome. Where the model disagrees with a published measurement, the model is
presumed wrong. Full references: REFERENCES.md.

## Mixed synapses annotated as chemical in EM
MaleCNS v1.0 lists 679 direct JO-A/B -> DNp01 chemical synapses. Physiology says this
contact is a mixed synapse that is primarily electrical (Pezier & Blagburn 2013).
`MIXED_IN_EM` declares such pairs; when the electrical model is on, the EM chemical
edges for those pairs are removed so the contact is not counted twice. The auditory
benchmark reports the EM-as-annotated arm for reference.

## MaleCNS naming notes (discovered on v1.0)
- No sugar/bitter GRN split exists. The taste afferents feeding MN9's excitatory inputs
  are discovered at runtime (sensory-superclass types with >= 50 synapses onto them;
  on v1.0 this is dominated by `GNG642`, with `BM_Taste` minor). Bitter suppression is SKIP.
- Shiu's FlyWire second-order names (Fdg, Bract, Zorro...) do not exist; the benchmark
  uses MN9's strongest input types (GNG###, DNge###) discovered at runtime and lists them.
- Sensory somas (JO, taste) have no `somaSide`; side is inferred from soma x vs midline.

## Known placeholders
- Loom tuning weights are ordinal from Turner et al. 2022 Fig 3A; replace with per-glomerulus dF/F.
- dF/F to firing-rate mapping unknown; PEAK_HZ anchored to GF physiology.
- Modulatory amines (DA, OA, 5-HT) are sign 0 and not modelled.
