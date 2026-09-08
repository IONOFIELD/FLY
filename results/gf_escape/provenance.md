# Provenance: gf_escape

## Neurotransmitter sign map
{'acetylcholine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1, 'dopamine': 0, 'octopamine': 0, 'serotonin': 0}

## Electrical synapses added (invisible to EM)
- DNp01 -> TTMn (ipsilateral, spikelet 9.0): shakB gap junction, 1:1 relay; Tanouye & Wyman 1980; Allen et al. 2006
- DNp01 -> PSI (ipsilateral, spikelet 9.0): shakB gap junction, 1:1 relay; Allen et al. 2006; Phelan et al. 2008
- JO-A* -> DNp01 (ipsilateral, spikelet {'compound_mV': 3.0, 'at_hz': 150}): JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009
- JO-B* -> DNp01 (ipsilateral, spikelet {'compound_mV': 3.0, 'at_hz': 150}): JON->GF electrical; Pezier & Blagburn 2013; Yorozu 2009

## Intrinsic overrides
- DNp01.b_adapt_mV = 30.0: all-or-none GF response; von Reyn 2014; Ache 2019
- superclass:cb_motor.b_adapt_mV = 10.0: sustained MN rates < 100 Hz; Azevedo 2020; McKellar 2020
- superclass:vnc_motor.b_adapt_mV = 10.0: sustained MN rates < 100 Hz; Azevedo 2020

## Regional gain: SEZ x2.0 (fitted, benchmarks/fit_regional.py), other x1

## Loom tuning (PLACEHOLDER ordinal, Turner et al. 2022 Fig 3A)
{'LC4': 1.0, 'LPLC2': 1.0, 'LC6': 1.0, 'LC16': 0.8, 'LC26': 0.6, 'LPLC1': 0.6, 'LC9': 0.5, 'LC17': 0.5, 'LC12': 0.4}

## Checks
- gf_spikes_per_response_1_to_2: PASS
- gf_response_prob_0.5_to_1: PASS
- gf_latency_10_60ms: PASS
- ttmn_one_to_one: PASS
- ttmn_lag_0.5_1.5ms: PASS
- cns_quiet: PASS
- electrical_necessary: PASS
- null_silent: PASS