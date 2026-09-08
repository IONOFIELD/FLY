# flycns: benchmarked spiking simulation of the Drosophila central nervous system

One-command, repeatable simulation of the Janelia FlyEM male CNS connectome
(v1.0, 166,700 neurons, brain + ventral nerve cord) as a signed leaky
integrate-and-fire network, validated against published circuit physiology.

## Launch

    ./fly

That is the whole interface: it activates the environment, loads your neuprint token from
`~/.neuprint_token`, and opens the menu. To run it from anywhere: `alias flycns=~/flycns/fly` in
your shell profile.

First time, from a fresh clone (macOS or Linux, Python 3.10 or newer):

    git clone https://github.com/IONOFIELD/flycns.git && cd flycns
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    echo "PASTE_YOUR_NEUPRINT_TOKEN" > ~/.neuprint_token      # from neuprint.janelia.org, account page
    ./fly

`./fly` opens the menu. On first run choose **1** (fetch connectome, ~5 min) and **2** (fetch
sensory annotations), then **4** (run the full suite, ~25 min). After that:

    16   3D cascade with real neuron morphology (opens in your browser)
    12   the same cascade drawn in the terminal
    11   interactive: choose a sensory group, stimulate it, watch what fires
    14   one neuron's skeleton and synapses in the terminal
     6   suite summary: every check, one parameter set

Every later session: `./fly` (or just `flycns` with the alias).

Without the menu, the same things are single commands: `./run_all.sh 20`, `python viz/cascade_3d.py`,
`python viz/explore.py --sez 2.0`, `python viz/neuron_ascii.py DNp01 --synapses`.

## Layout
    flycns/graph.py      load, NT signing, electrical synapses, intrinsic overrides, null model
    flycns/model.py      Brian2 LIF network (chemical + electrical synapses, adaptation)
    flycns/protocols.py  stimulus protocols with citations (loom after Turner, Krieger et al. 2022)
    benchmarks/          one file per circuit; writes results/<name>/report.json + provenance.md
    viz/cascade_3d.py    animated 3D cascade from saved spikes (plotly, html)
    viz/explore.py       INTERACTIVE: pick a sensory group and rate, drive it, watch the cascade
                         (one build, then ~10 s per stimulation; --sez GAIN to use a fitted SEZ gain)
    viz/neuron_ascii.py  one neuron's skeleton in braille (dorsal/side/front), with input and
                         output synapses and synapses onto its top partner; skeletons cached
    viz/cascade_ascii.py braille dorsal-view animation of the cascade in the terminal;
                         --synapses draws arbors and lights transmission sites (run
                         viz/fetch_cascade_synapses.py first, pulls synapse xyz from neuprint)
    make_testdata.py     synthetic graph for smoke tests (no data required)
    export_brainsets.py  simulated loom/sound/taste sessions as brainsets-style HDF5 for POYO+
                         (units table carries connectome features to seed unit embeddings)

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
searched class-wise gains (sensory, relay, local): taste -> MN9 never activated even at 3x
sensory gain, and the local gain feeding needs breaks the escape circuit. `benchmarks/fit_regional.py`
scales only synapses from SEZ intrinsic types (GNG, SAD, PRW, FLA, CAN) and selects gustatory
afferents by MaleCNS `subclass` (run `fetch_annotations.py` first).

`benchmarks/summarize.py` writes `results/SUITE.md`: every check across every circuit under one declared parameter set.

## Declared parameter set (v0.4)
Shiu 2024 LIF parameters; chemical kick 0.275 x 0.3 mV; SEZ intrinsic synapses x2.0 (fitted,
`benchmarks/fit_regional.py`); GF adaptation 30 mV (all-or-none, von Reyn 2014); motor neuron
adaptation 10 mV (sustained MN rates < 100 Hz, Azevedo 2020, McKellar 2020); electrical
synapses as listed. Every benchmark runs under this one set; `results/SUITE.md` reports it.

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

## Robustness (results/overnight, 2026-09-08)
Full suite under seeds 0-2, 5% and 10% neurotransmitter sign inversion, and w_scale 0.25-0.40.
At the fitted parameters, 18-19/19 across seeds; the only seed-sensitive checks were two
marginal specifications (GF response probability at exactly 0.80; summation tested on a
loom that already fired GF), both since re-specified with citations and dated notes in the
benchmark headers. 5% sign inversion costs 0-2 checks, 10% costs 4. Above w_scale 0.3 the
failures are flood signatures (extra TTMn spikes, indirect JO->GF depolarisation, extra-SEZ
activity); below it the suite holds to 0.25.

## Tier 2: fitting to recordings
The Turner, Krieger, Pang & Clandinin 2022 glomerulus imaging data are public (Dryad
doi:10.5061/dryad.h44j0zpp8). Menu 21 downloads them, menu 22 extracts the measured per-glomerulus
loom amplitudes and trial-gain distribution (`fit/extract_turner_loom.py`) and, when activated,
`flycns/protocols.py` uses them instead of the ordinal placeholder; provenance records the switch.

## Known placeholders
- Loom tuning weights are ordinal until `fit/extract_turner_loom.py` has been run on the Dryad data.
- dF/F to firing-rate mapping unknown; PEAK_HZ anchored to GF physiology.
- Modulatory amines (DA, OA, 5-HT) are sign 0 and not modelled.
