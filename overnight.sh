#!/usr/bin/env bash
# Robustness battery: seed replicates, neurotransmitter sign-flip perturbation, weight-scale sweep.
# Each run = full three-benchmark suite (~25 min). Saves results/overnight/<tag>/SUITE.md and commits.
# Usage:  caffeinate -i ./overnight.sh 2>&1 | tee overnight.log
set -u
source .venv/bin/activate 2>/dev/null || true
mkdir -p results/overnight
run () {   # run TAG SEED SIGNFLIP WSCALE
  local tag=$1; export FLYCNS_SEED=$2 FLYCNS_SIGNFLIP=$3 FLYCNS_WSCALE=$4 FLYCNS_RUN_TAG=$tag
  echo "=== $(date '+%H:%M') start $tag (seed $2, signflip $3, wscale $4)"
  for b in gf_escape auditory feeding; do python benchmarks/$b.py 40 data > "results/overnight/${tag}_${b}.log" 2>&1 || echo "  $b errored (see log)"; done
  python benchmarks/summarize.py > /dev/null 2>&1
  mkdir -p "results/overnight/$tag" && cp results/SUITE.md results/*/report.json "results/overnight/$tag/" 2>/dev/null
  grep -m1 "checks pass" results/SUITE.md
  git add -A results/overnight >/dev/null 2>&1 && git commit -qm "overnight: $tag $(grep -m1 -o '[0-9]*/[0-9]* checks pass' results/SUITE.md)" >/dev/null 2>&1 && git push -q 2>/dev/null
  echo "=== $(date '+%H:%M') done  $tag"
}
run seed1        1 0    0.3
run seed2        2 0    0.3
run signflip05_a 0 0.05 0.3
run signflip05_b 1 0.05 0.3
run signflip10   0 0.10 0.3
run wscale025    0 0    0.25
run wscale035    0 0    0.35
run wscale040    0 0    0.40
export FLYCNS_SEED=0 FLYCNS_SIGNFLIP=0 FLYCNS_WSCALE=0.3 FLYCNS_RUN_TAG=""
echo; echo "SUMMARY"; for d in results/overnight/*/; do printf "%-16s %s\n" "$(basename $d)" "$(grep -m1 -o '[0-9]*/[0-9]* checks pass' $d/SUITE.md 2>/dev/null)"; done
