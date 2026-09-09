"""
Neuropil ROI membership and meshes from the MaleCNS ROI volumes.

Two uses:
  1. a principled regional partition for the regional synaptic gain. The current
     SEZ definition is a type-name prefix heuristic (GNG/SAD/PRW/FLA/CAN); this
     replaces it with actual neuropil membership, using each neuron's synapse
     distribution across compartments (neuprint roiInfo, or the ROI volume).
  2. neuropil surface meshes for viz/cascade_3d.py, instead of the soma cloud.

Run:  python fit/roi_membership.py            (uses neuprint roiInfo; needs NEUPRINT_TOKEN)
      python fit/roi_membership.py --meshes   (also builds meshes from the ROI volumes;
                                               needs: pip install cloud-volume trimesh)
Writes data/roi_membership.parquet (bodyId, dominant ROI, fraction of synapses in it)
and, with --meshes, data/meshes/<roi>.npz for the 3D viewer.
"""
import os, sys, json
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("data"); OUT.mkdir(exist_ok=True)
SEZ_ROIS = ["GNG", "PRW", "SAD", "FLA", "CAN", "AMMC"]   # subesophageal zone compartments

if "--meshes" not in sys.argv:
    from neuprint import Client, fetch_custom
    tok = os.environ.get("NEUPRINT_TOKEN")
    if not tok:
        raise SystemExit('export NEUPRINT_TOKEN=$(cat ~/.neuprint_token) first')
    c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=tok)
    print("fetching per-neuron ROI synapse distributions...")
    q = "MATCH (n:Neuron) RETURN n.bodyId AS bodyId, n.roiInfo AS roi"
    df = fetch_custom(q, client=c)
    print(f"  {len(df):,} rows; roiInfo dtype example: {type(df.roi.iloc[0]).__name__}")
    rows, skipped = [], 0
    for bid, ri in zip(df.bodyId, df.roi):
        if isinstance(ri, str):                      # neuprint returns roiInfo as a JSON string
            try:
                ri = json.loads(ri)
            except Exception:
                skipped += 1; continue
        if not isinstance(ri, dict):
            skipped += 1; continue
        tot = {k: (v.get("pre", 0) + v.get("post", 0)) for k, v in ri.items() if isinstance(v, dict)}
        tot = {k: v for k, v in tot.items() if v > 0 and k not in ("CentralBrain", "VNC", "OpticLobe")}
        if not tot:
            continue
        s = sum(tot.values()); dom = max(tot, key=tot.get)
        sez = sum(v for k, v in tot.items() if any(k.startswith(r) for r in SEZ_ROIS))
        rows.append(dict(bodyId=bid, dominant_roi=dom, dominant_frac=tot[dom] / s,
                         sez_frac=sez / s, n_syn=s))
    if not rows:
        print(f"  no usable roiInfo ({skipped:,} skipped). First value was:\n  {str(df.roi.iloc[0])[:300]}")
        raise SystemExit("adjust the roiInfo parser to that shape")
    mem = pd.DataFrame(rows)
    print(f"  parsed {len(mem):,} neurons ({skipped:,} without usable roiInfo)")
    mem.to_parquet(OUT / "roi_membership.parquet", index=False)
    n_sez = int((mem.sez_frac > 0.5).sum())
    print(f"wrote data/roi_membership.parquet: {len(mem):,} neurons, {n_sez:,} with >50% of synapses in SEZ compartments")
    print(mem.dominant_roi.value_counts().head(15).to_string())
    print("\nflycns/graph.py can now define the SEZ region by membership (sez_frac > 0.5) "
          "instead of the type-prefix heuristic; compare the two sets before switching.")
else:
    try:
        from cloudvolume import CloudVolume
        import trimesh
    except ImportError:
        raise SystemExit("pip install cloud-volume trimesh")
    from skimage import measure
    MD = OUT / "meshes"; MD.mkdir(exist_ok=True)
    for name, path in [("brain", "gs://flyem-male-cns/rois/fullbrain-roi-v4"),
                       ("vnc", "gs://flyem-male-cns/rois/malecns-vnc-neuropil-roi-v0")]:
        # ROI volumes are sparse: absent chunks are normal, and full resolution is far
        # more than a surface needs, so read a coarse mip with fill_missing.
        v = CloudVolume(f"precomputed://{path}", use_https=True, progress=True,
                        fill_missing=True, mip=-1)
        print(name, "volume shape at coarsest mip", v.shape, "resolution", v.resolution)
        vol = np.asarray(v[:, :, :, 0]).astype(np.uint8)
        occupied = (vol > 0).astype(np.uint8)
        if occupied.sum() == 0:
            print("  volume empty at this mip; skipping"); continue
        verts, faces, _, _ = measure.marching_cubes(occupied, level=0.5, step_size=2)
        res = np.array(v.resolution)
        verts = verts * res
        np.savez_compressed(MD / f"{name}.npz", v=verts.astype("f4"), f=faces.astype("i4"))
        print(f"  wrote data/meshes/{name}.npz: {len(verts):,} vertices")
    print("viz/cascade_3d.py will use these if present")
