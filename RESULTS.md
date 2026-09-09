# Results note: MaleCNS v1.0 as a signed leaky integrate-and-fire network

*Status v0.6, 9 September 2026. A working note, not a manuscript. Everything here is
reproducible from `./run_all.sh` and `./overnight.sh`; numbers are from `results/`.*

## What was done
The complete male Drosophila central nervous system connectome (Janelia FlyEM MaleCNS v1.0,
166,700 neurons, brain and ventral nerve cord, released June 2026) was simulated as a
network of leaky integrate-and-fire neurons using the model and parameters of Shiu et al.
2024 (Nature), which validated the same approach on the female FlyWire brain. Synaptic
signs come from the connectome's neurotransmitter predictions. Three circuits with
published physiology were turned into pass/fail tests with null and ablation arms, and the
model was held to all three under one parameter set.

## What passes (18/18 scored checks, seeds 0-2)
**Loom escape.** Driving the loom-sensitive lobula columnar populations (tuning after Turner,
Krieger, Pang & Clandinin 2022) fires the giant fiber 1 to 2 times per responding trial
(1.06 mean) with ~30 ms latency and 0.80 response probability (95% CI 0.58-0.92). Each GF
spike drives exactly one TTMn spike 0.9 ms later through a modelled gap junction; without
the electrical model TTMn is silent; on a degree-preserving rewired connectome nothing fires.
This is a sensory to motor cascade across the neck connective at physiological latency.

**Auditory input to the giant fiber.** Johnston's organ drive alone leaves GF subthreshold
(~1 mV net), as recorded in vivo. Removing both the electrical model and the 679 EM synapses
annotated as chemical on this contact abolishes the response, identifying those synapses
as the pathway.

**Feeding: the model does not reproduce it, and we can say why.** Under the escape circuit's
parameters MN9 stays silent, and four uniform manipulations fail for one reason:

- MN9's synaptic input is balanced to 1.4% (2,966 excitatory vs 3,049 inhibitory synapses), so
  scaling anything moves both sides together.
- Gustatory afferents contact MN9's inhibitory relays three times more strongly than its
  excitatory ones (1,097 vs 3,333 synapses for the screen-selected subset; 1,314 vs 4,659 for
  all anatomically gustatory afferents).
- Signed path products from those afferents to MN9 are negative at two hops and positive at
  three to five: **the pathway is disinhibitory.**
- A silent network cannot express disinhibition. A uniform tonic depolarisation that makes it
  non-silent lets MN9 respond, but the wind-sensitive control drives it just as strongly
  (specificity lost) and the escape response is abolished.

Reproducing feeding therefore needs cell-specific spontaneous activity, that is, the
state-dependent modulation this model omits; Shiu et al. 2024 reported the same limitation for
inhibitory neurons in their 0 Hz-baseline model. Cell-type-specific spontaneous firing is
measured, not hypothetical: slow leg motor neurons fire at tens of Hz at rest, maintained by
cholinergic drive (blocking nicotinic receptors lowered both the rate and resting muscle force),
while fast motor neurons are silent (Azevedo et al. 2020, Fig. 3D and Fig. 4C). That is exactly
the pattern a uniform baseline cannot produce and is the form the missing mechanism should take. Full numbers and the scripts that produced
them: `results/feeding/mechanism.json`. What the feeding benchmark still scores: each of MN9's
excitatory second-order inputs drives it when stimulated directly, MN9 is silent on a rewired
null, activity stays inside the SEZ, and motor rates stay physiological.

## Auditory input to GF: contact set corrected from the connectome's own annotation
MaleCNS annotates JO neurons by modality (`subclass`: auditory / wind_gravity / grooming), and
that annotation does not follow the JO-A/JO-B type split: JO-B2, B3 and B4 are mostly
wind_gravity. The only JO types with chemical contacts onto DNp01 are JO-B1_a (541 synapses,
13 cells) and JO-B1_c (138, 6), both auditory; JO-A contributes none. The declared JON-GF
electrical synapse and the benchmark's stimulus sets now follow the modality annotation rather
than type-name prefixes, which previously mixed ~30 wind_gravity cells into the "sound" drive.
Kamikouchi et al. 2009 place the GF dendrite in AMMC zone A; since MaleCNS's A/B labels are not
the modality split, the two are not directly comparable and this is recorded as a naming
difference between datasets rather than a contradiction.

## Multisensory interaction at GF (reported, currently unresolved)
Whether auditory drive raises or lowers GF's response to a loom depends entirely on where the
loom sits relative to threshold, and earlier runs reported both directions at a fixed loom gain
(suppression 0.55 to 0.20 under the ordinal tuning; facilitation 0.00 to 0.07 under the measured
tuning at the same gain). Neither is quotable. The auditory benchmark now calibrates a
near-threshold working point per run (largest loom gain with GF hit rate <= 0.3 and mean
depolarisation >= 2 mV) and reports the effect there with a confidence interval on the
difference. Across the eight-run robustness battery at 20 trials per arm the reading was
facilitation four times, suppression twice and no change twice, all differences of 1 to 5
trials: the effect is NOT resolved at that trial count and no direction should be quoted.
Resolving a 0.15 to 0.30 difference needs roughly 120 trials per arm
(`FLYCNS_A2_TRIALS=120 python benchmarks/auditory.py`). The in vivo direction is not
established in our references either.

## Locomotor state changes the loom input (9 Sept 2026)
The datafiles carry the authors' own 50 Hz walking classification, so trials can be split by
behavioural state without touching the video archive. Labelling a loom trial walking if the fly
moves for more than half the stimulus window and stationary if less than a tenth (8 and 13 series
qualify), every glomerulus except LC12 is relatively larger when the fly walks: LC6 0.64 vs 0.40,
LC21 0.59 vs 0.39, LC26 0.93 vs 0.69, LC4 0.73 vs 0.52, LPLC2 1.00 vs 0.83. Trial-gain sigma is
also higher when walking (0.44 vs 0.36). This is the locomotor enhancement Turner, Krieger, Pang
& Clandinin 2022 report, recovered by an independent extraction. Because the benchmarks simulate
a fly that is not walking, the **stationary** amplitudes are the default input;
`FLYCNS_LOOM_STATE=walking` or `=all` selects the others.

## Measured loom input (added 8 Sept 2026)
The ordinal loom tuning has been replaced by amplitudes extracted from the public Turner,
Krieger, Pang & Clandinin 2022 dataset (Dryad doi:10.5061/dryad.h44j0zpp8): 10 flies, 150 loom
trials, 13 glomeruli, peak dF/F normalised to the strongest glomerulus (LC17 1.00, LC12 0.75,
LC26 0.74, LPLC2 0.69, LPLC1 0.51, LC4 0.50, LC16 0.49, LC6 0.40; small-object types 0.29-0.37),
and a measured trial-gain sigma of 0.38 (assumed 0.5 before). Under measured input, GF response
probability falls from 0.80 to 0.50, still within von Reyn's range. The dF/F-to-rate scale
(5 Hz for the strongest type) is the one remaining free parameter of the stimulus and is
declared as such.

## A general limitation: this model class cannot express disinhibition
Three independent circuits in this connectome turn out to work by sign inversion, and a network
with no spontaneous activity cannot represent any of them:
- **Taste to MN9**: signed path products negative at two hops, positive at three to five.
- **Lamina to Mi1 (the canonical ON pathway)**: Mi1's largest input is L1 with 141,873
  glutamatergic (inhibitory) synapses. The ON response is a double inversion, light depolarising
  photoreceptors, histamine inhibiting L1, L1 releasing less glutamate and Mi1 being released.
- **AstA release**: Pm3 is GABAergic, so activating it in a silent network produces nothing at
  all, in the model or in principle.
This is why the feeding failure is not a quirk of one pathway. Any circuit whose output is
carried by the removal of inhibition is invisible to a zero-baseline LIF, and a uniform baseline
does not fix it (it destroys stimulus specificity and abolishes escape;
`benchmarks/fit_baseline.py`). The missing ingredient is cell-specific spontaneous activity,
which is measured for at least one motor pool (slow leg MNs fire at tens of Hz at rest,
cholinergically maintained; Azevedo et al. 2020).

## Peptidergic modulation, declared (AstA/Pm3/AstA-R1)
The feeding result concludes that the missing mechanism is cell-specific modulation. One instance
is fully specified in the literature and visible in this connectome: Krieger 2023 (PhD thesis,
Stanford, Ch. 2) shows AstA is released by a single visual cell type, Pm3, that AstA-R1 is
expressed by L1-L5, C2, Mi1, Mi15, Dm9, Tm2, TmY3 and T2, that AstA perfusion or optogenetic Pm3
activation increases the peak-to-trough dynamic range of Mi1's light response, and that AstA-R1
knockdown in Mi1 blocks it while leaving conventional transmission intact. All 12 receptor-
expressing types and 69 Pm3 cells are present in MaleCNS v1.0, where Pm3 also inhibits Mi1 with
18,394 GABAergic synapses: the two channels run between the same cells, which is why the biology
needed perfusion and knockdown to separate them. `flycns/graph.py` declares this as
`PEPTIDE_MODULATION` (source, receptor-expressing targets, gain, timescale, citation), off by
default, and `benchmarks/peptide_asta.py` runs the model mirror of the genetic experiment
(synaptic only / plus AstA / AstA with Mi1 excluded). The result: **not testable in this model class.** Mi1 fires
0.0006 spikes per cell per flash under 100 Hz drive to its largest cholinergic inputs, so a
multiplicative gain has nothing to act on and all three modulated arms are identical. Mi1 is
disinhibition-driven (363,509 inhibitory against 132,350 excitatory synapses; its largest input
is inhibitory L1), and Pm3 is GABAergic, so activating it in a silent network can do nothing
even in principle. The layer is kept as declared, cited infrastructure, and the arms are not
tuned to make Mi1 fire. Beyond that, the measured effect is a change in a graded biphasic calcium
waveform, and every cell in this model spikes, so even a working version could only offer a
spike-count analogue.

**One blocking item now accounts for three separate results.** Feeding, the ON pathway and AstA
modulation all fail because the same ingredient is missing: cell-specific spontaneous activity.
That makes it the next thing to build, ahead of any fourth circuit.

## What the model needed that the connectome does not contain
- Electrical synapses (invisible to EM): GF to TTMn, GF to PSI, JON to GF. Documented
  anatomy, added with citations and calibrated to measured potentials.
- (Retired) a regional gain on subesophageal-zone synapses. Fitted to 2.0 while the SEZ was
  defined by type-name prefixes, where it appeared to open the feeding pathway. Defining the
  SEZ from the connectome's own compartment annotations (>50% of synapses in
  GNG/PRW/SAD/FLA/CAN/AMMC/PENP, excluding sensory/motor/efferent: 3,519 neurons, 4.0% of
  synaptic weight) shows the population is net inhibitory onto that pathway under either
  definition, and MN9 stays silent at every gain. The parameter is now 1.0 and the apparent
  effect is documented as an artefact of the naming heuristic. Without
  it, taste never reaches the motor neurons at gains where the escape circuit is stable;
  with a global gain, escape and feeding cannot both pass. Stated as a hypothesis about
  SEZ synapse strength, fitted, and open to test.
- Spike-frequency adaptation on GF and motor neurons. The GF increment (30 mV) is the
  smallest value satisfying von Reyn's constraints in a pre-specified sweep
  (`benchmarks/fit_gf_adapt.py`); without it GF fires up to 8 spikes per strong loom.
  Motor neuron adaptation (proboscis pool only) keeps sustained rates under 100 Hz (Azevedo
  2020, McKellar 2020); relay motor neurons (TTMn) are not adapted, as they follow GF 1:1
  (Tanouye & Wyman 1980).
- A rescale of Shiu's unitary weight to 0.3x, needed once the nerve cord and its hub
  neurons are included.

## What could not be tested
MaleCNS v1.0 does not annotate taste modality. Sugar versus bitter, and therefore Shiu's
contralateral MN9 prediction for labellar sugar GRNs and any bitter-suppression test, are
recorded as not applicable rather than passed or failed. Driving every gustatory sensillum
at once yields no feeding output, consistent with a mixed sugar and bitter population.

## How much rests on uncertain transmitter calls
MaleCNS publishes per-T-bar transmitter probabilities alongside the aggregate per-neuron label
the model signs with. Over the annotated neurons that carry the graph, mean per-synapse
confidence is 0.78; neurons whose individual T-bars mostly disagree with their aggregate label
carry **3.5% of total synaptic weight**, and low-confidence neurons (mean < 0.6) carry 0.15%.
`FLYCNS_SIGNFLIP_TARGETED=1` inverts exactly those neurons, which is a sharper test than the
random 5%/10% flips and is included in the robustness battery (`fit/synapse_nt.py`).

## Robustness
Escape passes 8/8 unchanged when every neuron whose per-T-bar transmitter predictions disagree
with its aggregate label is sign-inverted (132,680 edges, 2.0% of synaptic weight): the circuit
does not rest on any uncertain transmitter call.

Benchmarks now default to 40 trials: at 20, GF response probability (criterion 0.5) failed on
the two runs where it sat at the floor with a confidence interval spanning it, which is a
resolution limit rather than a model result. Full suite at seeds 0-2: 18/18. Inverting 5% of neurotransmitter signs at random costs
0-2 checks; 10% costs 4. The shared weight scale holds from 0.25 to 0.30 and floods above.

## Does the cord inherit the wrong parameters? (bounded, not resolved)
Shiu et al.'s single-neuron values were chosen for central-brain neurons and we apply them to
the nerve cord too. `benchmarks/sensitivity_cord.py` sweeps cord-only tau_m (5-40 ms), threshold
gap (4-10 mV) and refractory period (1-5 ms) across all 14,153 cord neurons. **Every escape
check holds at every setting, and every metric is invariant** (GF 1.06 spikes per response,
response probability 0.68, TTMn/GF ratio 1.06, relay lag 0.90 ms, population rate 0.0045 Hz;
standard deviation 0.0 across the sweep). The escape result therefore does not depend on cord
single-neuron parameters: the cascade crosses the neck through a declared electrical relay and
one strong chemical connection, neither of which the membrane time constant gates. The inherited
brain parameters remain a limitation for any circuit that asks the cord to compute, and are not
one for what this benchmark tests. Azevedo et al. 2020 (eLife 9:e56754) measured the relevant properties for tibia flexor motor
neurons: input resistance 150 MOhm (fast), 300 (intermediate), 700 (slow), with resting
potential, spontaneous rate and soma, neurite and axon diameter covarying along the same
gradient. Those numbers are recorded in `flycns/graph.py` but are NOT entered as parameters,
for three stated reasons: tau_m = R*C and capacitance is not reported, so the time constant is
not derivable; using the gradient as a postsynaptic gain requires a fast/slow assignment for
MaleCNS motor neurons that does not exist (the 2024 FANC connectome identified leg MNs by
muscle target in a different dataset); and somatic current injection failed to evoke spikes in
fast and intermediate MNs because the spike initiation zone is electrically isolated from the
soma, which a single-compartment cell cannot represent at all.

## What this does not claim
The single-neuron parameters (membrane time constant, threshold, refractory period, unitary
weight) are Shiu et al.'s values for central-brain neurons of a female fly, applied uniformly
to every neuron of a male brain and nerve cord. Neurotransmitter and cell-type labels are
MaleCNS's own, not transferred. Sexual dimorphism is unlikely to matter for these circuits
(the Cell paper reports 8,069 isomorphic types; none of the circuit neurons here is a known
dimorphic type), but the uniform-membrane assumption is coarser in the cord, where motor
neurons are large, low-resistance cells; the motor-neuron adaptation override partly
compensates for that. Per-superclass membrane parameters constrained by motor neuron
recordings (Azevedo 2020, McKellar 2020) are the natural next refinement.

No parameter here was fitted to a neural recording; "passes" means consistent with published
measurements under the stated assumptions. The regional gain is a fitted number with a
plausible story. The loom tuning is ordinal pending the authors' per-glomerulus data.
Modulatory transmitters are omitted. Where the model disagrees with a measurement, the
model is presumed wrong.

## Relation to prior work
Model class, parameters, and feeding test design: Shiu et al. 2024. Escape physiology:
Tanouye & Wyman 1980, Allen et al. 2006, von Reyn et al. 2014/2017, Ache et al. 2019.
Auditory-GF synapse: Pezier & Blagburn 2013, Yorozu et al. 2009. Loom tuning: Turner et al.
2022. What is added: the test harness, the multi-circuit single-parameter-set constraint,
the brain-plus-cord spiking cascade on MaleCNS, and the documented negatives. Full list in
REFERENCES.md.
