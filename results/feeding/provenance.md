# Provenance: feeding (port of Shiu et al. 2024 tests to MaleCNS)

sugar types used: ['aPhM2a', 'aPhM5', 'PhG1c', 'claw_tpGRN', 'LB3a', 'LB3b', 'LB3c', 'LB3d']
bitter types used: none resolved
second-order present: ['DNge062', 'GNG120', 'GNG117', 'GNG108']
motor present: ['MN9', 'MN6', 'MN8']

## Checks
- F2b_extra_SEZ_quiet: PASS
- F2c_motor_rate_physiological: PASS
- F2_monotonic: PASS
- F1_contra_gt_ipsi: FAIL
- F3_bitter_suppresses: SKIP
- F4_second_order_sufficient: PASS
- F5_null_silent: PASS