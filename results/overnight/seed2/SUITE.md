# MaleCNS v1.0 LIF benchmark suite  [seed2]

env: FLYCNS_SEED=2 FLYCNS_SIGNFLIP=0 FLYCNS_WSCALE=0.3

**18/19 checks pass under ONE parameter set** (same LIFParams for every circuit).

## Parameters
```
{
  "v_rest_mV": -52.0,
  "v_thresh_mV": -45.0,
  "v_reset_mV": -52.0,
  "tau_m_ms": 20.0,
  "tau_syn_ms": 5.0,
  "refractory_ms": 2.2,
  "w_syn_mV": 0.275,
  "w_scale": 0.3,
  "tau_adapt_ms": 100.0,
  "gap_delay_ms": 0.8,
  "chem_delay_ms": 1.8,
  "dt_ms": 0.1,
  "seed": 2,
  "class_gains": {
    "sensory": 1.0,
    "relay": 1.0,
    "local": 1.0
  },
  "region_gains": {
    "SEZ": 2.0,
    "other": 1.0
  }
}
```
effective chemical kick per synapse: 0.0825 mV (Shiu 2024 unitary 0.275 mV x MaleCNS rescale 0.3)

## Checks
| benchmark | check | result |
|---|---|---|
| auditory | A1_jo_alone_mostly_subthreshold | PASS |
| auditory | A2_summation_with_loom | FAIL |
| auditory | A3_declared_pairs_carry_jo_drive | PASS |
| auditory | A4_null_silent | PASS |
| auditory | A5_cns_quiet | PASS |
| feeding | F2b_extra_SEZ_quiet | PASS |
| feeding | F2c_motor_rate_physiological | PASS |
| feeding | F2_monotonic | PASS |
| feeding | F1_laterality_matches_wiring | PASS |
| feeding | F1_shiu_contralateral_bias | SKIP |
| feeding | F3_bitter_suppresses | SKIP |
| feeding | F4_second_order_sufficient | PASS |
| feeding | F5_null_silent | PASS |
| gf_escape | gf_spikes_per_response_1_to_2 | PASS |
| gf_escape | gf_hit_ge_0.8 | PASS |
| gf_escape | gf_latency_10_60ms | PASS |
| gf_escape | ttmn_one_to_one | PASS |
| gf_escape | ttmn_lag_0.5_1.5ms | PASS |
| gf_escape | cns_quiet | PASS |
| gf_escape | electrical_necessary | PASS |
| gf_escape | null_silent | PASS |

Declared deviations from the raw connectome: see each results/*/provenance.md and flycns/graph.py.

## Citation check: 25/25 tags resolve in REFERENCES.md