# References

Every literature tag used in `flycns/` and `benchmarks/` resolves here. Entries marked
[verify] have a DOI or details that should be confirmed against the source before citing
in a manuscript. This project builds on the work below; it does not contradict it.

## Connectome and simulation framework
- **MaleCNS v1.0** — Janelia FlyEM Project Team et al. A connectome of the male Drosophila central nervous system. *Cell* (2026), published 3 Sept 2026. Data: male-cns.janelia.org, neuprint `male-cns:v1.0`, CC-BY 4.0. [verify: final author list and DOI]
- **Shiu et al. 2024** — Shiu PK, Sterne GR, Spiller N, Blagburn JM, et al. A Drosophila computational brain model reveals sensorimotor processing. *Nature* 634:210–219 (2024). doi:10.1038/s41586-024-07763-9. Source of LIF parameters (V_rest/V_reset −52 mV, V_th −45 mV, tau_m 20 ms, tau_syn 5 ms, refractory 2.2 ms, delay 1.8 ms, 0.275 mV/synapse), the feeding and grooming validation design, and the treatment of neurotransmitter sign.
- **Shiu et al. 2022** — Shiu PK, Sterne GR, Engert S, Dickson BJ, Scott K. Taste quality and hunger interactions in a feeding sensorimotor circuit. *eLife* 11:e79887 (2022). Feeding circuit cell types (Fdg, Bract, Zorro, Roundup, Usnea, Phantom, Clavicle, Rattle, Cleaver).
- **Lappalainen et al. 2024** — Lappalainen JK, Tschopp FD, Prakhya S, et al. Connectome-constrained networks predict neural activity across the fly visual system. *Nature* 634:1132–1140 (2024). doi:10.1038/s41586-024-07939-3. Template for differentiable parameter fitting. [verify DOI]
- **Eckstein et al. 2024** — Eckstein N, Bates AS, Champion A, et al. Neurotransmitter classification from electron microscopy images at synaptic sites in Drosophila melanogaster. *Cell* 187:2574–2594 (2024). Basis of predicted neurotransmitter labels.
- **Brian2** — Stimberg M, Brette R, Goodman DFM. Brian 2, an intuitive and efficient neural simulator. *eLife* 8:e47314 (2019).

## Neurotransmitter sign
- **Liu & Wilson 2013** — Liu WW, Wilson RI. Glutamate is an inhibitory neurotransmitter in the Drosophila olfactory system. *PNAS* 110:10294–10299 (2013).
- **Hardie 1989** — Hardie RC. A histamine-activated chloride channel involved in neurotransmission at a photoreceptor synapse. *Nature* 339:704–706 (1989).

## Giant fiber escape circuit
- **Tanouye & Wyman 1980** — Tanouye MA, Wyman RJ. Motor outputs of giant nerve fiber in Drosophila. *J Neurophysiol* 44:405–421 (1980). GF→TTM ~0.8–0.9 ms, GF→DLM ~1.4 ms latencies.
- **Allen et al. 2006** — Allen MJ, Godenschwege TA, Tanouye MA, Phelan P. Making an escape: development and function of the Drosophila giant fibre system. *Semin Cell Dev Biol* 17:31–41 (2006). GF–TTMn and GF–PSI shakB gap junctions.
- **Phelan et al. 2008** — Phelan P, Goulding LA, Tam JLY, et al. Molecular mechanism of rectification at identified electrical synapses in the Drosophila giant fiber system. *Curr Biol* 18:1955–1960 (2008).
- **von Reyn et al. 2014** — von Reyn CR, Breads P, Peek MY, et al. A spike-timing mechanism for action selection. *Nat Neurosci* 17:962–970 (2014). GF fires at most a single spike per loom; most looms leave GF subthreshold; single spike suffices for takeoff.
- **von Reyn et al. 2017** — von Reyn CR, Nern A, Williamson WR, et al. Feature integration drives probabilistic behavior in the Drosophila escape response. *Neuron* 94:1190–1204 (2017). LC4 and LPLC2 inputs to GF; ~39–42° size threshold.
- **Ache et al. 2019** — Ache JM, Polsky J, Alghailani S, et al. Neural basis for looming size and velocity encoding in the Drosophila giant fiber escape pathway. *Curr Biol* 29:1073–1081 (2019). 55 LC4 + 108 LPLC2 synapse onto GF.
- **Augustin et al. 2019** — Augustin H, Zylbertal A, Partridge L. A computational model of the escape response latency in the giant fiber system of Drosophila melanogaster. *eNeuro* 6(2):ENEURO.0423-18.2019. Fitted gap-junction conductance model.

## Auditory / mechanosensory input to GF
- **Kamikouchi et al. 2009** — Kamikouchi A, Inagaki HK, Effertz T, et al. The neural basis of Drosophila gravity-sensing and hearing. *Nature* 458:165–171 (2009). JO-A/B vibration-sensitive subgroups.
- **Yorozu et al. 2009** — Yorozu S, Wong A, Fischer BJ, et al. Distinct sensory representations of wind and near-field sound in the Drosophila brain. *Nature* 458:201–205 (2009).
- **Pezier & Blagburn 2013 (Pézier)** — Pézier A, Blagburn JM. Auditory responses of engrailed and empty spiracles mutants in Drosophila; JON→GF transmission is predominantly electrical, abolished in shakB2. *PLoS ONE* (2013). [verify: exact title and volume]

## Auditory pathway physiology
- **Tootoonian et al. 2012** — Tootoonian S, Coen P, Kawadler JM, Murthy M. Neural representations of courtship song in the Drosophila brain. *J Neurosci* 32:787-798 (2012). AMMC-A1, B1, B2 respond to sound with graded, non-spiking potentials; GF spikes to somatic current injection and sound summates with subthreshold injection to produce a full spike (Fig. 2D); GF is not activated by visual or mechanosensory stimuli alone; AMMC-B1 cholinergic, B2 GABAergic; AMMC-A1 and GF gap-junction coupled.
- **Vaughan et al. 2014** — Vaughan AG, Zhou C, Manoli DS, Baker BS. Neural pathways for the detection and discrimination of conspecific song in Drosophila melanogaster. *Curr Biol* 24:1039-1049 (2014). aPN1 (= AMMC-B1) and aLN(al) are the only AMMC cell types necessary for behavioural courtship-song responses in either sex; aPN1 integrates pulse rate over 25-65 ms IPI with attenuation at 20 ms and an AMMC-to-WED band-pass near 35 ms; aPN1 reported non-spiking.
- **Kim et al. 2020** — Kim H, et al. Synaptic organisation of JO-A and JO-B axons (FAFB EM). *J Comp Neurol* (2020). JO-A axons longer with more presynaptic sites (896 +/- 344) than JO-B (329 +/- 89); both contact AMMC-B1 and GF; type-1 JO-A tuned to 100-200 Hz.

## Loom stimulus and LC neuron physiology
- **Turner, Krieger, Pang & Clandinin 2022** — Turner MH, Krieger A, Pang MM, Clandinin TR. Visual and motor signatures of locomotion dynamically shape a population code for feature detection in Drosophila. *eLife* 11:e82587 (2022). Pan-glomerulus imaging; loom-responsive LC types; shared trial gain.
- **Krieger 2023** — Krieger A. How do neuropeptides shape sensory processing in Drosophila. PhD thesis, Stanford University (2023). Ch. 2: AstA is released by a single visual-system cell type, Pm3; AstA-R1 expressed in L1-L5, C2, Mi1, Mi15, Dm9, Tm2, TmY3, T2 (AstA-R2 in L2); AstA perfusion and optogenetic Pm3 activation increase the peak-to-trough dynamic range of Mi1's response to light flashes; AstA-R1 knockdown in Mi1 blocks the increase and reduces the hyperpolarising phase, with conventional transmission intact.

## Nerve cord connectome and motor neuron identity
- **Azevedo et al. 2024** — Azevedo A, Lesser E, Phelps JS, Mark B, et al. Connectomic reconstruction of a female Drosophila ventral nerve cord. *Nature* (2024). Leg and wing motor neuron identification with muscle targets (69-70 MNs per T1 leg; SETi 7,090 and FETi 14,904 input synapses). Anatomy and identity, no biophysical measurements.

## Motor neuron physiology
- **Azevedo et al. 2020** — Azevedo AW, Dickinson ES, Gurung P, Cherry L, Ahmed AW, Tuthill JC. A size principle for recruitment of Drosophila leg motor neurons. *eLife* 9:e56754 (2020). Tibia flexor MN input resistance 150/300/700 MOhm for fast/intermediate/slow (Fig. 3E, n=15/11/14); covarying resting potential, spontaneous firing rate and soma, neurite and axon diameter (Fig. 3C-E); slow MN spontaneous firing maintained by cholinergic input and setting resting muscle force (Fig. 4C); somatic current injection fails to evoke spikes in fast and intermediate MNs (spike initiation zone electrically isolated from the soma); recruitment ordered weakest to strongest (Fig. 5).

## Feeding circuit
- **Gordon & Scott 2009** — Gordon MD, Scott K. Motor control in a Drosophila taste circuit. *Neuron* 61:373–384 (2009). MN9 necessary and sufficient for proboscis extension.
- **McKellar et al. 2020** — McKellar CE, Siwanowicz I, Dickson BJ, Simpson JH. Controlling motor neurons of every muscle for fly proboscis reaching. *eLife* 9:e54978 (2020).

## Larval connectome (referenced in planning)
- **Winding et al. 2023** — Winding M, Pedigo BD, Barnes CL, et al. The connectome of an insect brain. *Science* 379:eadd9330 (2023).

## Decoding (planned)
- **Azabou et al. 2023** — Azabou M, Arora V, Ganesh V, et al. A unified, scalable framework for neural population decoding (POYO). *NeurIPS* (2023). arXiv:2310.16046.
- **Azabou et al. 2025** — Azabou M, et al. Multi-session, multi-task neural decoding from distinct cell types and brain regions (POYO+). *ICLR* (2025).

## Prior MaleCNS dynamics demonstrations (acknowledged, not built on)
- Georgia Tech / E. Smith MaleCNS-in-Minecraft demonstration (Sept 2026, social media; code not public at time of writing).
- Eon Systems embodied FlyWire model (company writeup, March 2026).
