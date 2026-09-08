# Results note: MaleCNS v1.0 as a signed leaky integrate-and-fire network

*Status v0.4.2, 8 September 2026. A working note, not a manuscript. Everything here is
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

**Feeding.** Pharyngeal and labellar gustatory afferents drive the proboscis motor pool
(MN9, MN10, MN11, MN4, MN6, MN7; 11 of 43 central-brain motor types) monotonically with
input rate, at motor rates of ~70 Hz, with the rest of the CNS silent and a wind-sensitive
control population producing nothing. Unilateral drive produces ipsilateral MN9 output,
matching the 73:1 structural bias in the wiring.

## Multisensory interaction at GF (reported, direction withdrawn)
Whether auditory drive raises or lowers GF's response to a loom depends entirely on where the
loom sits relative to threshold, and earlier runs reported both directions at a fixed loom gain
(suppression 0.55 to 0.20 under the ordinal tuning; facilitation 0.00 to 0.07 under the measured
tuning at the same gain). Neither is quotable. The auditory benchmark now calibrates a
near-threshold working point per run (largest loom gain with GF hit rate <= 0.3 and mean
depolarisation >= 2 mV) and reports the effect there; the value is in
`results/auditory/report.json` under `A2_prediction` with its calibration table. The in vivo
direction is not established in our references.

## Measured loom input (added 8 Sept 2026)
The ordinal loom tuning has been replaced by amplitudes extracted from the public Turner,
Krieger, Pang & Clandinin 2022 dataset (Dryad doi:10.5061/dryad.h44j0zpp8): 10 flies, 150 loom
trials, 13 glomeruli, peak dF/F normalised to the strongest glomerulus (LC17 1.00, LC12 0.75,
LC26 0.74, LPLC2 0.69, LPLC1 0.51, LC4 0.50, LC16 0.49, LC6 0.40; small-object types 0.29-0.37),
and a measured trial-gain sigma of 0.38 (assumed 0.5 before). Under measured input, GF response
probability falls from 0.80 to 0.50, still within von Reyn's range. The dF/F-to-rate scale
(5 Hz for the strongest type) is the one remaining free parameter of the stimulus and is
declared as such.

## What the model needed that the connectome does not contain
- Electrical synapses (invisible to EM): GF to TTMn, GF to PSI, JON to GF. Documented
  anatomy, added with citations and calibrated to measured potentials.
- A regional gain: synapses from subesophageal-zone intrinsic neurons scaled x2.0. Without
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

## Robustness
Full suite at seeds 0-2: 18/18. Inverting 5% of neurotransmitter signs at random costs
0-2 checks; 10% costs 4. The shared weight scale holds from 0.25 to 0.30 and floods above.

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
