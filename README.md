# flycns: benchmarked spiking simulation of the Drosophila central nervous system

A leaky integrate-and-fire network built from the complete male fly connectome (Janelia FlyEM
MaleCNS v1.0: 166,700 neurons, brain and ventral nerve cord, released June 2026), held to
published circuit physiology through a test suite with pass/fail criteria, null models and
ablation arms. Three circuits pass under one declared parameter set. Every deviation from the
raw wiring is cited. Status: v0.4.3, September 2026. Results note: `RESULTS.md`.

## Launch

    ./flycns

That is the whole interface: it activates the environment, loads your neuprint token from
`~/.neuprint_token`, and opens the menu. To run it from anywhere: `alias flycns=~/flycns/fly`.

First time, from a fresh clone (macOS or Linux, Python 3.10 or newer):

    git clone https://github.com/IONOFIELD/FLYCNS.git flycns && cd flycns
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    echo "PASTE_YOUR_NEUPRINT_TOKEN" > ~/.neuprint_token      # neuprint.janelia.org, account page
    ./fly

In the menu: **1** fetch the connectome (~5 min), **2** fetch sensory annotations, **4** run the
suite (~25 min). Then **16** for the 3D cascade with real morphology, **11** to stimulate a sensory
group and watch what fires, **14** for one neuron's skeleton and synapses, **6** for the summary.

## What passes (18/18 scored checks, seeds 0 to 2)

| circuit | stimulus | readout | criteria (published physiology) |
|---|---|---|---|
| loom escape | 9 loom-sensitive LC types (Turner et al. 2022) | DNp01 (GF), TTMn, PSI | 1-2 GF spikes per response, P(response) 0.5-1.0 with CI, 10-60 ms latency, 1:1 TTMn relay at 0.8 ms, relay lost without electrical synapse, rewired null silent, CNS quiet |
| auditory to GF | JO-A/JO-B at 150 Hz | DNp01 membrane V | sound alone subthreshold, pathway carried by the declared JON-GF contact, null silent, CNS quiet |
| feeding | pharyngeal/labellar gustatory afferents | proboscis motor neurons | monotonic dose response, MN rates < 100 Hz, activity confined to SEZ, unilateral drive matches 73:1 structural laterality, second-order sufficiency, null silent |

Recorded but not scored: **the model predicts auditory drive suppresses GF response to a
near-threshold loom (0.55 to 0.20-0.40)**; not established in vivo in our references. Not
applicable with this annotation: bitter suppression and Shiu's labellar-sugar contralateral bias
(MaleCNS does not annotate taste modality).

Robustness (`results/overnight`): 5% neurotransmitter sign inversion costs 0-2 checks, 10% costs 4;
the shared weight scale holds from 0.25 to 0.30 and floods above.

## Declared parameter set

Shiu et al. 2024 LIF parameters; chemical kick 0.275 mV x 0.3 (MaleCNS rescale, fitted by sweep);
SEZ intrinsic synapses x2.0 (fitted, `benchmarks/fit_regional.py`); GF adaptation 30 mV (smallest
value satisfying von Reyn 2014 constraints, `benchmarks/fit_gf_adapt.py`); motor neuron adaptation
10 mV (sustained rates < 100 Hz); electrical synapses GF-TTMn, GF-PSI (1:1 relay, 0.8 ms) and
JON-GF (calibrated to a 3 mV compound potential). The 679 EM synapses annotated as chemical on the
JON-GF contact are removed when the electrical model is on, so the contact is not counted twice.
Monoamines are sign 0 by deliberate scope. Full list with citations: `flycns/graph.py`, `REFERENCES.md`.

## Fitting to recordings (in progress)

The Turner, Krieger, Pang & Clandinin 2022 glomerulus imaging data are public (Dryad
doi:10.5061/dryad.h44j0zpp8; Dryad blocks scripted downloads, so menu **21** gives browser
instructions). Menu **22** runs `fit/extract_turner_loom.py`, which produces measured per-glomerulus
loom amplitudes and the trial-gain distribution; when activated, `flycns/protocols.py` uses them in
place of the ordinal placeholder and provenance records the switch. Next: per-type LC gain fit so
the simulated population reproduces the measured relative amplitudes.

## Layout

    fly / fly.py          launcher and menu
    flycns/graph.py       load, NT signing, electrical synapses, intrinsic and regional overrides, null model
    flycns/model.py       Brian2 LIF (chemical + electrical synapses, adaptation, class and regional gains)
    flycns/protocols.py   stimulus protocols with provenance (loom after Turner et al. 2022)
    flycns/anatomy.py     braille projection of the CNS for terminal views
    benchmarks/           gf_escape, auditory, feeding; fit_gains, fit_regional, fit_gf_adapt;
                          screen_afferents, diagnose_feeding; summarize -> results/SUITE.md
    fit/                  extraction of measured tuning from published recordings
    viz/                  cascade_3d (skeletons, synapses, camera presets), cascade_ascii, explore
                          (interactive), neuron_ascii (single cell), fetch_cascade_synapses
    export_brainsets.py   simulated sessions as brainsets-style HDF5 for POYO+, with connectome
                          features per unit
    overnight.sh          robustness battery: seeds, sign flips, weight scale
    RESULTS.md            results note; REFERENCES.md full references

## MaleCNS naming notes (discovered on v1.0)
- No sugar/bitter GRN split. Gustatory afferents are anatomical subclasses: `labellar bristle`,
  `taste peg`, `pharyngeal sensillum`. Driving all at once gives no feeding output (expected for a
  mixed population); the scored arm uses the screen-selected MN9-driving subset, labelled as such.
- Shiu's FlyWire feeding names (Fdg, Bract, ...) do not exist; MN9's inputs are `GNG###`, `DNge###`.
- Sensory somas lie outside the CNS; side comes from the `_L/_R` instance suffix.

## What this builds on and what it adds
A port and extension of Shiu et al. 2024 (Nature) from FlyWire to MaleCNS. The model, its
parameters and the feeding validation design are theirs; escape physiology from von Reyn, Ache,
Tanouye & Wyman, Allen, Augustin; the JON-GF synapse from Pezier & Blagburn; loom tuning from
Turner, Krieger, Pang & Clandinin. Added here: the test harness with numeric criteria and nulls,
one parameter set across independent circuits, a spiking sensory-to-motor cascade across the neck
connective on MaleCNS, the documented negatives, and one testable prediction. Where the model
disagrees with a measurement, the model is presumed wrong.

## Known limits
Single-neuron parameters are Shiu's central-brain values applied uniformly to brain and cord; no
parameter is yet fitted to a recording (the fitting work above changes that); modulatory transmitters omitted;
gap junctions only where documented.
