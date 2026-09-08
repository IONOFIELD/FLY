#!/usr/bin/env python
"""
flycns launcher: one menu for setup, benchmarks, fits, exploration and visualization.

    python fly.py
"""
import json, os, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
STATE = ROOT / ".flycns_state.json"
state = json.loads(STATE.read_text()) if STATE.exists() else {"cell": "DNp01"}

def save():
    STATE.write_text(json.dumps(state))

def run(*cmd, env=None):
    e = dict(os.environ); e.update(env or {})
    try:
        subprocess.run([sys.executable, *cmd] if cmd[0].endswith(".py") else list(cmd), env=e)
    except KeyboardInterrupt:
        print("\n(interrupted)")

def status():
    have = lambda p: "yes" if Path(p).exists() else "no "
    t = os.environ.get("NEUPRINT_TOKEN", "")
    tok = "set" if (len(t) > 40 and "<" not in t and "token" not in t.lower()) else \
          "NOT SET or placeholder. Save it once:  echo PASTE > ~/.neuprint_token ; then  export NEUPRINT_TOKEN=$(cat ~/.neuprint_token)"
    suite = ""
    if Path("results/SUITE.md").exists():
        for line in Path("results/SUITE.md").read_text().splitlines():
            if "checks pass" in line:
                suite = line.strip("* "); break
    print(f"""
  data/neurons.parquet {have('data/neurons.parquet')}   data/edges.parquet {have('data/edges.parquet')}   annotations {have('data/neurons.parquet') if 'subclass' in _cols() else 'no '}
  neuprint token: {tok}
  last suite: {suite or 'none run'}
  selected cell for neuron view: {state['cell']}""")

def _cols():
    try:
        import pyarrow.parquet as pq
        return pq.read_schema("data/neurons.parquet").names
    except Exception:
        return []

def ask(prompt, default=""):
    v = input(f"{prompt}{' [' + default + ']' if default else ''}> ").strip()
    return v or default

MENU = """
  flycns  ::  benchmarked spiking simulation of the Drosophila CNS (MaleCNS v1.0)
  ------------------------------------------------------------------------------
  SETUP          1  fetch connectome (neurons + edges, ~5 min)
                 2  fetch sensory annotations (subclass, entry nerve)
                 3  status
  BENCHMARKS     4  run full suite (3 circuits, ~25 min)      5  run one benchmark
                 6  show suite summary                        7  robustness battery (overnight)
  FITS           8  GF adaptation sweep       9  SEZ regional gain sweep      10  afferent screen
  EXPLORE       11  interactive: pick a sensory group, stimulate, watch the cascade
  VISUALIZE     12  whole-CNS cascade (braille, terminal)     13  same, with synapse sites
                14  single neuron view  (current: {cell})     15  choose cell for neuron view
                16  3D cascade (html, opens in browser)
  EXPORT        20  simulated sessions in brainsets/POYO+ layout (HDF5 with connectome unit features)
  DOCS          17  results note      18  references      19  provenance of last run
                 q  quit
"""

while True:
    print(MENU.format(cell=state["cell"]))
    c = ask("select").lower()
    if c in ("q", "quit", ""):
        break
    elif c == "1": run("fetch_malecns.py")
    elif c == "2": run("fetch_annotations.py")
    elif c == "3": status()
    elif c == "4": run("bash", "run_all.sh", ask("trials", "20"))
    elif c == "5":
        b = ask("benchmark (gf_escape / auditory / feeding)", "gf_escape")
        run(f"benchmarks/{b}.py", ask("trials", "20"))
        if b == "gf_escape": state["cell"] = "DNp01"; save()
        if b == "feeding": state["cell"] = "MN9"; save()
    elif c == "6":
        run("benchmarks/summarize.py")
    elif c == "7":
        print("  runs ~2.5 h; use:  caffeinate -i ./overnight.sh 2>&1 | tee overnight.log")
        if ask("start now? (y/n)", "n").lower() == "y": run("bash", "overnight.sh")
    elif c == "8": run("benchmarks/fit_gf_adapt.py", ask("trials", "20"))
    elif c == "9": run("benchmarks/fit_regional.py", ask("trials", "6"))
    elif c == "10":
        run("benchmarks/screen_afferents.py", ask("trials", "3"), "data", ask("rate Hz", "50"), ask("w_scale", "0.3"))
    elif c == "11":
        run("viz/explore.py", "--sez", ask("SEZ gain", "2.0"))
    elif c == "12":
        run("viz/cascade_ascii.py", *([] if ask("layout: landscape / portrait", "landscape") == "landscape" else ["--portrait"]))
    elif c == "13":
        if not Path("results/gf_escape/cascade_sites.parquet").exists():
            print("  fetching synapse sites for the cascade neurons (needs token)..."); run("viz/fetch_cascade_synapses.py")
        run("viz/cascade_ascii.py", "--synapses")
    elif c == "14":
        flags = []
        v = ask("view: dorsal / side / front", "dorsal")
        if v != "dorsal": flags.append(f"--{v}")
        if ask("show synapses? (y/n)", "y").lower() == "y": flags.append("--synapses")
        if ask("mark synapses onto top partner? (y/n)", "n").lower() == "y": flags.append("--partners")
        run("viz/neuron_ascii.py", state["cell"], *flags)
    elif c == "15":
        state["cell"] = ask("cell type or bodyId (e.g. DNp01, TTMn, MN9, LC4, 10001)", state["cell"]); save()
    elif c == "16":
        run("viz/cascade_3d.py"); run("open", "results/gf_escape/cascade_3d.html")
    elif c == "20": run("export_brainsets.py", ask("trials per session", "40"))
    elif c == "17": run("less", "RESULTS.md")
    elif c == "18": run("less", "REFERENCES.md")
    elif c == "19":
        for p in Path("results").glob("*/provenance.md"): print(f"\n--- {p}\n" + p.read_text())
    else:
        print("  unknown option")
    input("\n  enter to return to menu...")
