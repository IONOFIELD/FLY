# Provenance: feeding (port of Shiu et al. 2024 tests to MaleCNS)

sugar types used: ['LB1a', 'LB1b', 'LB1c', 'LB1d', 'LB1e', 'LB2a', 'LB2b', 'LB2c', 'LB2d', 'LB3', 'LB3a', 'LB3b', 'LB3c', 'LB3d', 'LB4a', 'LB4b', 'PhG10', 'PhG11', 'PhG12', 'PhG13', 'PhG14', 'PhG15', 'PhG16', 'PhG1a', 'PhG1b', 'PhG1c', 'PhG2', 'PhG3', 'PhG4', 'PhG5', 'PhG6', 'PhG7', 'PhG8', 'PhG9', 'aPhM1', 'aPhM2a', 'aPhM2b', 'aPhM3', 'aPhM4', 'aPhM5', 'claw_tpGRN', 'dorsal_tpGRN']
bitter types used: none resolved
second-order present: ['DNge062', 'GNG120', 'GNG117', 'GNG108']
motor present: ['MN9', 'MN6', 'MN8']

## Checks
- F1_contra_gt_ipsi: FAIL
- F2b_extra_SEZ_quiet: PASS
- F2c_motor_rate_physiological: PASS
- F2_monotonic: FAIL
- F3_bitter_suppresses: SKIP
- F4_second_order_sufficient: PASS
- F5_null_silent: PASS