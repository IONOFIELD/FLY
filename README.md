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

`benchmarks/summarize.py` writes `results/SUITE.md`: every check across every circuit under one declared parameter set.

## Every deviation from the raw connectome is declared
See `flycns/graph.py` (SIGN_MAP, ELECTRICAL_SYNAPSES, INTRINSIC_OVERRIDES) and
`results/*/provenance.md`. Effective chemical kick is 0.275 mV (Shiu 2024 unitary
weight) x 0.3 (MaleCNS rescale from sweep) = 0.0825 mV per synapse; chemical delay
1.8 ms (Shiu 2024). Electrical synapses: 1:1 relays use a 9 mV spikelet with 0.8 ms
latency (GF-TTMn); population inputs (JON-GF) are calibrated to a measured compound
potential and share one coupling budget so resting partners cannot shunt the target.
Monoamines are sign 0 by deliberate scope (Shiu used +1; no principled default exists).

## What is prior art and what is new
Shiu et al. 2024 validated feeding and grooming on FlyWire with fixed literature
parameters; the feeding benchmark here is a port. The new claims are: a benchmark
suite with numeric pass/fail criteria, one parameter set across independent circuits,
and a validated spiking sensory-to-motor cascade crossing the neck connective on MaleCNS.

## MaleCNS naming notes (discovered on v1.0)
- Labellar taste GRNs are one type, `BM_Taste`; there is no sugar/bitter split, so the
  bitter-suppression check is skipped and the feeding benchmark reports it as SKIP.
- Shiu's FlyWire second-order names (Fdg, Bract, Zorro...) do not exist; the benchmark
  uses MN9's strongest input types (GNG###, DNge###) discovered at runtime and lists them.
- Sensory somas (JO, taste) have no `somaSide`; side is inferred from soma x vs midline.

## Known placeholders
- Loom tuning weights are ordinal from Turner et al. 2022 Fig 3A; replace with per-glomerulus dF/F.
- dF/F to firing-rate mapping unknown; PEAK_HZ anchored to GF physiology.
- Modulatory amines (DA, OA, 5-HT) are sign 0 and not modelled.
