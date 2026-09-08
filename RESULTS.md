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

## What the model predicts (not scored)
Auditory drive during a near-threshold loom **reduces** GF response probability from 0.55 to
0.20-0.40 (three seeds). The electrical JON to GF path is excitatory but small; chemical
auditory pathways recruit inhibition onto GF that outweighs it. We found no in vivo
measurement of this interaction in our references. It is a testable prediction.

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
  Motor neuron adaptation keeps sustained rates under 100 Hz (Azevedo 2020, McKellar 2020).
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
