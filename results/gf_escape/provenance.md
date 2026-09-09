# Provenance: gf_escape

## Neurotransmitter sign map
{'acetylcholine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1, 'dopamine': 0, 'octopamine': 0, 'serotonin': 0}

## Electrical synapses added (invisible to EM)
- DNp01 -> TTMn (ipsilateral, spikelet 9.0): shakB gap junction, 1:1 relay; Tanouye & Wyman 1980; Allen et al. 2006
- DNp01 -> PSI (ipsilateral, spikelet 9.0): shakB gap junction, 1:1 relay; Allen et al. 2006; Phelan et al. 2008
- JO-A* -> DNp01 (ipsilateral, spikelet {'compound_mV': 3.0, 'at_hz': 150}): JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009
- JO-B* -> DNp01 (ipsilateral, spikelet {'compound_mV': 3.0, 'at_hz': 150}): JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009

## Intrinsic overrides
- DNp01.b_adapt_mV = 15.0: all-or-none GF response; von Reyn 2014; Ache 2019. Fitted by pre-stated rule (smallest value giving 1-2 spikes/response, <=2 on strong looms, P(response) 0.5-1, TTMn relay intact; benchmarks/fit_gf_adapt.py). 30 mV under the ordinal loom tuning (2026-09-08), 15 mV under the measured Turner 2022 tuning (2026-09-08); 0 mV gives up to 8 spikes.
- superclass:cb_motor.b_adapt_mV = 10.0: sustained MN rates < 100 Hz; Azevedo 2020; McKellar 2020

## Regional gain: SEZ x1.0 (fitted, benchmarks/fit_regional.py), other x1
   SEZ definition: anatomical (roiInfo, 3,519 SEZ-intrinsic neurons)

## Loom tuning: MEASURED (all trials): Turner, Krieger, Pang & Clandinin 2022, eLife 11:e82587; Dryad doi:10.5061/dryad.h44j0zpp8 (21 flies, extracted 2026-09-09); trial gain sigma 0.3687; peak rate 5.0 Hz (free dF/F->rate scale)
{'LC17': 1.0, 'LC12': 0.9416, 'LPLC2': 0.8695, 'LC26': 0.8051, 'LC15': 0.6731, 'LC4': 0.5948, 'LPLC1': 0.5639, 'LC16': 0.537, 'LC6': 0.448, 'LC21': 0.4218, 'LC18': 0.4135, 'LC9': 0.3532, 'LC11': 0.307}

## Checks
- gf_spikes_per_response_1_to_2: PASS
- gf_response_prob_0.5_to_1: PASS
- gf_latency_10_60ms: PASS
- ttmn_one_to_one: PASS
- ttmn_lag_0.5_1.5ms: PASS
- cns_quiet: PASS
- electrical_necessary: PASS
- null_silent: PASS