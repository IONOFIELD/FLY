"""
Fetch bulk MaleCNS files from Google Storage (no token needed; CC-BY).

  synapse-nt   per-T-bar neurotransmitter probabilities (2.7 GB) -> per-synapse confidence
               and sign, instead of one aggregate label per presynaptic neuron
  roi-brain    brain neuropil compartment volume (256 nm) -> neuropil membership + meshes
  roi-vnc      VNC neuropil compartment volume (256 nm)

Run:  python fetch_bulk.py synapse-nt          (or roi-brain, roi-vnc, all)
Files land in data/bulk/. Resumable.
"""
import subprocess, sys
from pathlib import Path

BASE = "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome"
TARGETS = {
    "synapse-nt": (f"{BASE}/tbar-neurotransmitters-male-cns-v1.0.feather", "tbar-neurotransmitters.feather", "2.7 GB"),
    "annotations": (f"{BASE}/body-annotations-male-cns-v1.0-minconf-0.5.feather", "body-annotations.feather", "13 MB"),
    "connectome": (f"{BASE}/connectome-weights-male-cns-v1.0-minconf-0.5.feather", "connectome-weights.feather", "1.1 GB"),
}
OUT = Path("data/bulk"); OUT.mkdir(parents=True, exist_ok=True)
which = sys.argv[1] if len(sys.argv) > 1 else "synapse-nt"
keys = list(TARGETS) if which == "all" else [which]
for k in keys:
    if k not in TARGETS:
        raise SystemExit(f"unknown target {k}; choose from {list(TARGETS)} or 'all'")
    url, name, size = TARGETS[k]
    print(f"fetching {k} ({size}) -> data/bulk/{name}")
    subprocess.run(["curl", "-L", "-C", "-", "--progress-bar", "-o", str(OUT / name), url])
print("""
ROI volumes are neuroglancer-precomputed on gs:// and need cloud-volume:
    pip install cloud-volume
    python -c "from cloudvolume import CloudVolume; v=CloudVolume('precomputed://gs://flyem-male-cns/rois/fullbrain-roi-v4', use_https=True); print(v.shape)"
See fit/roi_membership.py once downloaded.""")
