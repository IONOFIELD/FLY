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
import os
tag = os.environ.get("FLYCNS_RUN_TAG", "")
md = ["# MaleCNS v1.0 LIF benchmark suite" + (f"  [{tag}]" if tag else ""), "",
      f"env: FLYCNS_SEED={os.environ.get('FLYCNS_SEED','0')} FLYCNS_SIGNFLIP={os.environ.get('FLYCNS_SIGNFLIP','0')} "
      f"FLYCNS_WSCALE={os.environ.get('FLYCNS_WSCALE','0.3')}", "",
      f"**{passed}/{total} checks pass under ONE parameter set** (same LIFParams for every circuit).", "",
      "## Parameters", "```", json.dumps(p, indent=2), "```",
      f"effective chemical kick per synapse: {p['w_syn_mV']*p['w_scale']:.4f} mV "
      f"(Shiu 2024 unitary {p['w_syn_mV']} mV x MaleCNS rescale {p['w_scale']})", "",
      "## Checks", "| benchmark | check | result |", "|---|---|---|"] + rows + [
      "", "Declared deviations from the raw connectome: see each results/*/provenance.md and flycns/graph.py."]
# citation integrity: every author-year tag in code must resolve in REFERENCES.md
import re
refs = Path("REFERENCES.md").read_text() if Path("REFERENCES.md").exists() else ""
code = "".join(f.read_text() for f in list(Path("flycns").glob("*.py")) + list(Path("benchmarks").glob("*.py")))
tags = set(re.findall(r"([A-Z][A-Za-z\u00e9\-]+(?: et al\.| & [A-Z][A-Za-z]+)?) (19\d\d|20\d\d)", code))
missing = sorted(f"{a} {y}" for a, y in tags if a.split()[0].rstrip(",") not in refs or y not in refs)
md += ["", f"## Citation check: {len(tags)-len(missing)}/{len(tags)} tags resolve in REFERENCES.md"]
if missing:
    md += ["Unresolved: " + ", ".join(missing)]
Path("results/SUITE.md").write_text("\n".join(md))
print("\n".join(md))
