"""Collect every results/*/report.json into results/SUITE.md with one parameter set declared."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import LIFParams
rows, total, passed = [], 0, 0
for rep in sorted(Path("results").glob("*/report.json")):
    r = json.loads(rep.read_text())
    for k, v in r.get("checks", {}).items():
        total += v is not None; passed += bool(v)
        rows.append(f"| {r['benchmark']} | {k} | {'SKIP' if v is None else ('PASS' if v else 'FAIL')} |")
p = LIFParams().__dict__
md = ["# MaleCNS v1.0 LIF benchmark suite", "",
      f"**{passed}/{total} checks pass under ONE parameter set** (same LIFParams for every circuit).", "",
      "## Parameters", "```", json.dumps(p, indent=2), "```",
      f"effective chemical kick per synapse: {p['w_syn_mV']*p['w_scale']:.4f} mV "
      f"(Shiu 2024 unitary {p['w_syn_mV']} mV x MaleCNS rescale {p['w_scale']})", "",
      "## Checks", "| benchmark | check | result |", "|---|---|---|"] + rows + [
      "", "Declared deviations from the raw connectome: see each results/*/provenance.md and flycns/graph.py."]
Path("results/SUITE.md").write_text("\n".join(md))
print("\n".join(md))
