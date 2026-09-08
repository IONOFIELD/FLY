"""
Pull synapse coordinates for the neurons active in a cascade window, so the
animation can show arbors and the actual transmission sites.

Reads  results/gf_escape/spikes_real.parquet (+ trial choice as cascade_ascii)
Writes results/gf_escape/cascade_sites.parquet     every pre/post site of active neurons (arbors)
       results/gf_escape/cascade_links.parquet     synapses between active pre -> active post

Run:  python viz/fetch_cascade_synapses.py [trial]      (needs NEUPRINT_TOKEN)
"""
import os, sys
import numpy as np
import pandas as pd
from pathlib import Path
from neuprint import Client, fetch_synapses, fetch_synapse_connections, NeuronCriteria as NC

token = os.environ.get("NEUPRINT_TOKEN")
if not token:
    raise SystemExit('Set NEUPRINT_TOKEN first: export NEUPRINT_TOKEN="..."')
c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=token)

RES = Path("results/gf_escape")
sp = pd.read_parquet(RES / "spikes_real.parquet")
trials = pd.read_csv(RES / "real_electrical.csv")
TRIAL = int(sys.argv[1]) if len(sys.argv) > 1 else int((trials.gf - 1.5).abs().idxmin())
onset = 300 + TRIAL * 360
win = sp[(sp.t_ms >= onset) & (sp.t_ms < onset + 120)]
ids = sorted(win.bodyId.unique().tolist())
# always include the escape circuit output cells so their arbors are drawn even if silent
neurons = pd.read_parquet("data/neurons.parquet").set_index("bodyId")
ids = sorted(set(ids) | set(neurons.index[neurons.type.isin(["DNp01", "TTMn", "PSI"])]))
print(f"trial {TRIAL}: {len(ids)} neurons; fetching synapse sites...")

sites = fetch_synapses(NC(bodyId=ids), client=c)
sites = sites[["bodyId", "type", "x", "y", "z"]].rename(columns={"type": "kind"})   # kind: pre | post
sites["celltype"] = sites.bodyId.map(neurons["type"]); sites["superclass"] = sites.bodyId.map(neurons["superclass"])
sites.to_parquet(RES / "cascade_sites.parquet", index=False)
print(f"  {len(sites):,} synaptic sites on {sites.bodyId.nunique()} neurons")

links = fetch_synapse_connections(NC(bodyId=ids), NC(bodyId=ids), client=c)
links = links[["bodyId_pre", "bodyId_post", "x_pre", "y_pre", "z_pre", "x_post", "y_post", "z_post"]]
links["pre_type"] = links.bodyId_pre.map(neurons["type"]); links["post_type"] = links.bodyId_post.map(neurons["type"])
links.to_parquet(RES / "cascade_links.parquet", index=False)
print(f"  {len(links):,} synapses between active neurons; top pairs:")
print(links.groupby(["pre_type", "post_type"]).size().sort_values(ascending=False).head(10).to_string())
