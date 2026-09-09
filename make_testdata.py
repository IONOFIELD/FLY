"""Tiny synthetic graph with the GF circuit types so the pipeline can be smoke-tested without MaleCNS."""
import numpy as np, pandas as pd
rng = np.random.default_rng(0)
rows, bid = [], 1000
def add(t, n, sc, side_alt=True):
    global bid
    for k in range(n):
        rows.append(dict(bodyId=bid, type=t, superclass=sc, somaSide=("L","R")[k%2] if side_alt else "L",
                         consensusNt="acetylcholine", x=rng.normal(), y=rng.normal(), z=rng.normal()))
        bid += 1
for t in ["LC4","LPLC2","LC6","LC16","LC26","LPLC1","LC9","LC17","LC12"]: add(t, 20, "visual_projection")
add("DNp01", 2, "descending_neuron"); add("TTMn", 2, "vnc_motor"); add("PSI", 2, "vnc_intrinsic")
add("JO-A1", 20, "cb_sensory"); add("JO-B1", 20, "cb_sensory")
add("Gr64f_sugar_GRN", 20, "cb_sensory"); add("Gr66a_bitter_GRN", 10, "cb_sensory")
add("GNG001", 4, "cb_intrinsic"); add("Bract", 4, "descending_neuron"); add("MN9", 2, "cb_motor")
add("JO-FV", 10, "cb_sensory"); add("aPhM2a", 5, "cb_sensory"); add("BM_Taste", 10, "cb_sensory")
add("Pm3", 6, "ol_intrinsic"); add("Mi1", 20, "ol_intrinsic"); add("Tm2", 10, "ol_intrinsic")
add("filler", 300, "cb_intrinsic")
n = pd.DataFrame(rows)
n["subclass"] = None
n.loc[n.type.str.startswith("JO-A")|n.type.str.startswith("JO-B"), "subclass"] = "auditory"
n.loc[n.type=="JO-FV", "subclass"] = "wind_gravity"
n.loc[n.type=="aPhM2a", "subclass"] = "pharyngeal sensillum"; n.loc[n.type=="Gr64f_sugar_GRN", "subclass"] = "taste peg"
n.loc[n.type.str.startswith("JO-") | n.type.str.contains("GRN"), "somaSide"] = None
n["instance"] = n.type + "_" + n.somaSide.fillna(pd.Series(["L","R"]*(len(n)//2+1))[:len(n)].set_axis(n.index))
n.loc[n.type.str.startswith("JO-"), "instance"] = None
e = []
gf = n[n.type=="DNp01"].bodyId.values
for b in n[n.type.isin(["LC4","LPLC2","LC6"])].bodyId: 
    for g in gf: e.append(dict(pre=b, post=g, weight=int(rng.integers(60,120)), nt="acetylcholine"))
def chain(pre_t, post_t, w, contra=False):
    for a in n[n.type==pre_t].itertuples():
        for b in n[n.type==post_t].itertuples():
            if contra and a.somaSide == b.somaSide: continue
            e.append(dict(pre=a.bodyId, post=b.bodyId, weight=w, nt="acetylcholine"))
chain("L1","Mi1",30); chain("Pm3","Mi1",20); chain("L2","Pm3",30); chain("Gr64f_sugar_GRN", "GNG001", 120); chain("aPhM2a", "GNG001", 150); chain("BM_Taste", "GNG001", 60); chain("GNG001", "MN9", 400, contra=True); chain("Bract", "MN9", 400)
for a in n[n.type=="Gr66a_bitter_GRN"].bodyId:
    for b in n[n.type=="Fdg"].bodyId: e.append(dict(pre=a, post=b, weight=40, nt="gaba"))
for a in n[n.type=="JO-B1"].bodyId[:5]:
    for g in gf: e.append(dict(pre=a, post=g, weight=100, nt="acetylcholine"))
for i in range(3000):
    a, b = rng.choice(n.bodyId.values, 2, replace=False)
    e.append(dict(pre=a, post=b, weight=int(rng.integers(5,10)), nt=rng.choice(["acetylcholine","gaba"])))
pd.DataFrame(e).to_parquet("testdata/edges.parquet", index=False); n.to_parquet("testdata/neurons.parquet", index=False)
print(len(n), "neurons", len(e), "edges")
