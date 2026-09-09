# flycns: benchmarked spiking simulation of the Drosophila central nervous system

A leaky integrate-and-fire network built from the complete male fly connectome (Janelia FlyEM
MaleCNS v1.0: 166,700 neurons, brain and ventral nerve cord, released June 2026), held to
published circuit physiology by a test suite with numeric pass/fail criteria, null models and
ablation arms. Three circuits, one declared parameter set, no fitted gains. Every deviation
from the raw wiring is cited. Status: v0.6, September 2026. Results note: `RESULTS.md`.

## Launch

    ./fly

That is the whole interface: it activates the environment, loads your neuprint token from
`~/.neuprint_token`, and opens the menu. To run it from anywhere: `alias flycns=~/flycns/fly`.

First time, from a fresh clone (macOS or Linux, Python 3.10 or newer):

    git clone https://github.com/IONOFIELD/FLYCNS.git flycns && cd flycns
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    echo "PASTE_YOUR_NEUPRINT_TOKEN" > ~/.neuprint_token      # neuprint.janelia.org, account page
    ./fly

In the menu: **1** fetch the connectome (~5 min), **2** fetch sensory annotations, **4** run the
suite (~40 min). Then **16** for the 3D cascade with real morphology, **11** to stimulate a
sensory group and watch what fires, **14** for one neuron's skeleton and synapses, **6** for the
summary, **27** for the measured account of why feeding fails.

## What the suite reports (16/16 scored, seeds 0 to 2)

| circuit | stimulus | readout | scored criteria |
|---|---|---|---|
| loom escape | 9 loom-sensitive LC types, amplitudes **measured** from Turner et al. 2022 | DNp01 (GF), TTMn, PSI | 1-2 GF spikes per response; response probability 0.5-1.0 with a Wilson CI; 10-60 ms latency; 1:1 TTMn relay at 0.8 ms; relay lost without the electrical synapse; rewired null silent; CNS quiet |
| auditory to GF | JO-A/JO-B at 150 Hz | DNp01 membrane potential | sound alone subthreshold; the drive is carried by the declared JON-GF contact; null silent; CNS quiet |
| feeding | gustatory afferents by anatomical subclass | proboscis motor neurons | each of MN9's excitatory second-order inputs drives it directly; MN9 silent on a rewired null; activity confined to the SEZ; motor rates < 100 Hz |

**Feeding does not reproduce sugar-driven extension, and the reason is measured.** MN9's input
is balanced to 1.4%; gustatory afferents contact its inhibitory relays about three times more
strongly than its excitatory ones; signed path products are negative at two hops and positive at
three to five, so the pathway is disinhibitory; a silent network cannot express that, and a
uniform tonic baseline that makes it non-silent destroys specificity and abolishes escape. Four
uniform manipulations were tested and all fail for this reason (`results/feeding/mechanism.json`,
menu **27**). Reproducing feeding needs cell-specific spontaneous activity, i.e. the
state-dependent modulation this model omits; Shiu et al. 2024 reported the same limitation.

Also reported, not scored: the effect of auditory drive on GF's response to a near-threshold
loom, at a working point calibrated per run. At 40 trials per arm the difference is
**statistically unresolved** (95% CI −0.33 to +0.23); `FLYCNS_A2_TRIALS=120` for a usable
estimate. Not applicable with this annotation: bitter suppression and Shiu's labellar-sugar
contralateral bias (MaleCNS does not annotate taste modality).

## Robustness

- Escape passes 8/8 **unchanged** when every neuron whose per-T-bar transmitter predictions
  disagree with its aggregate label is sign-inverted (132,680 edges, 2.0% of synaptic weight):
  the circuit does not rest on any uncertain transmitter call.
- Transmitter signing overall: mean per-synapse confidence 0.91 over the 173,308 annotated
  neurons; 1.2% of neurons below 0.6; neurons whose T-bars disagree with their label carry 3.5%
  of synaptic weight (`fit/synapse_nt.py`).
- Random sign inversion: 5% costs 0-2 checks, 10% costs 3-4. Seeds 0-2: 16-18/18 across versions.
- The shared weight scale holds from 0.25 to 0.30 and floods above (`overnight.sh`).

## Declared parameter set

Shiu et al. 2024 LIF parameters (tau_m 20 ms, threshold −45 mV, rest/reset −52 mV, refractory
2.2 ms, synaptic tau 5 ms, delay 1.8 ms); chemical kick 0.275 mV × 0.3 (MaleCNS rescale, fitted
by sweep); GF spike-triggered adaptation 15 mV (smallest value satisfying von Reyn 2014's
constraints under measured tuning, `benchmarks/fit_gf_adapt.py`); proboscis motor neuron
adaptation 10 mV (sustained rates < 100 Hz), relay motor neurons deliberately not adapted
(TTMn follows GF 1:1, Tanouye & Wyman 1980); electrical synapses GF-TTMn and GF-PSI (1:1 relay,
0.8 ms) and JON-GF (calibrated to a 3 mV compound potential); the 679 EM synapses annotated as
chemical on the JON-GF contact removed when the electrical model is on, so the contact is not
counted twice; monoamines sign 0 by deliberate scope; regional and class gains all 1.0; no tonic
baseline. Full list with citations: `flycns/graph.py`, `REFERENCES.md`.

## Measured inputs

The loom stimulus uses per-glomerulus amplitudes and a trial-gain distribution extracted from
the public Turner, Krieger, Pang & Clandinin 2022 dataset (Dryad doi:10.5061/dryad.h44j0zpp8;
10 flies, 150 loom trials): LC17 1.00, LC12 0.75, LC26 0.74, LPLC2 0.69, LPLC1 0.51, LC4 0.50,
LC16 0.49, LC6 0.40, small-object types 0.29-0.37, gain sigma 0.38. Menu **21** gives download
instructions (Dryad blocks scripted downloads), **22** runs the extraction
(`fit/extract_turner_loom.py`). The dF/F-to-rate scale (5 Hz for the strongest type) is the one
remaining free parameter of the stimulus and is declared as such.

## Layout

    fly / fly.py          launcher and menu
    flycns/graph.py       load, NT signing, electrical synapses, intrinsic and regional overrides,
                          null models, targeted sign perturbation
    flycns/model.py       Brian2 LIF (chemical + electrical synapses, adaptation, gains, optional baseline)
    flycns/protocols.py   stimulus protocols with provenance; measured loom tuning when available
    flycns/anatomy.py     braille projection of the CNS for terminal views
    benchmarks/           gf_escape, auditory, feeding; feeding_mechanism (the negative result);
                          fit_gains, fit_regional, fit_gf_adapt, fit_baseline; screen_afferents,
                          diagnose_feeding; summarize -> results/SUITE.md
    fit/                  extraction from published data: Turner loom tuning, per-synapse
                          transmitter confidence, neuropil ROI membership
    viz/                  cascade_3d (skeletons, synapses, camera presets), cascade_ascii,
                          explore (interactive), neuron_ascii, fetch_cascade_synapses
    export_brainsets.py   simulated sessions as brainsets-style HDF5 for POYO+, with connectome
                          features per unit
    overnight.sh          robustness battery: seeds, random and targeted sign flips, weight scale
    RESULTS.md            results note; REFERENCES.md full references

## MaleCNS naming notes (discovered on v1.0)
- No sugar/bitter GRN split. Gustatory afferents are anatomical subclasses: `labellar bristle`,
  `taste peg`, `pharyngeal sensillum`. Driving all at once gives no feeding output.
- Shiu's FlyWire feeding names (Fdg, Bract, ...) do not exist; MN9's inputs are `GNG###`, `DNge###`.
- Sensory somas lie outside the CNS; side comes from the `_L/_R` instance suffix.
- Region membership is available per neuron from `roiInfo` (`fit/roi_membership.py`), which is
  what retired an earlier type-name-prefix definition of the SEZ.

## What this builds on and what it adds
A port and extension of Shiu et al. 2024 (Nature) from FlyWire to MaleCNS. The model, its
parameters and the feeding validation design are theirs; escape physiology from von Reyn, Ache,
Tanouye & Wyman, Allen, Augustin; the JON-GF synapse from Pezier & Blagburn; loom tuning from
Turner, Krieger, Pang & Clandinin. Added here: the test harness with numeric criteria and nulls,
one parameter set across independent circuits, a spiking sensory-to-motor cascade across the neck
connective on MaleCNS, measured stimulus tuning from public imaging, an uncertainty-targeted
perturbation, and a mechanistic account of one circuit the model class cannot reproduce. Where
the model disagrees with a measurement, the model is presumed wrong.

## Known limits
Single-neuron parameters are Shiu's central-brain values applied uniformly to brain and cord;
no parameter is fitted to a neural recording (only the stimulus is measured); modulatory
transmitters are omitted, which is exactly what the feeding result implicates; gap junctions
appear only where documented.
