#!/usr/bin/env python3
"""Review existing projections, no fitting, no exclusions, no clinical releases."""
import csv,hashlib,math,sys,statistics
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
src=repo/"private_data/v8_joint461_fixed500_compatibility_v1_0_1"
prep=repo/"private_data/v8_joint461_full_G026_QC_normalization_v1"
fuscc=repo/"private_data/v8_fuscc_author_annotation_addendum_v1/360_fuscc_author_annotation_evidence_LOCAL_ONLY.tsv"
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def index(rows,key):
    d={r[key]:r for r in rows};assert len(d)==len(rows)
    return d
def num(x):
    v=float(x);assert math.isfinite(v);return v
def truth(x):
    s=x.strip().upper();assert s in {"TRUE","FALSE"};return s=="TRUE"
def stats(v):
    if not v:return "NA"
    v=sorted(v);i=(len(v)-1)*.95;lo=math.floor(i);hi=math.ceil(i)
    return {"median":round(statistics.median(v),7),"p95":round(v[lo]+(i-lo)*(v[hi]-v[lo]),7),"max":round(v[-1],7)}
def desc(rows):
    v=[r["tvd"] for r in rows]
    bad=[r for r in rows if not r["same"]]
    good=[r for r in rows if r["same"]]
    return {"n":len(v),"TVD":stats(v),"n_gt_0.05":sum(x>.05 for x in v),"n_gt_0.10":sum(x>.10 for x in v),
            "dominant_switches":len(bad),"median_CPM_margin_switch":round(statistics.median([r["margin"] for r in bad]),7) if bad else None,
            "median_CPM_margin_nonswitch":round(statistics.median([r["margin"] for r in good]),7) if good else None}
assert sha(fuscc)=="1cae682fe3ee006df02e22654f558c0ee2e28ab66fd5a422c67374609e1d1dd5"
F=index([r for r in read(fuscc) if r["membership"]=="RETAINED334"],"experiment_acc_LOCAL_ONLY")
assert len(F)==334
manifest=read(prep/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv")
assert len(manifest)==461 and all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in manifest)
expected={"SRP157974":334,"SRP057793":90,"SRP002628":14,"SRP094795":23}
assert Counter(r["study"] for r in manifest)==expected
for study,n in expected.items():
    ids={r["library_id_LOCAL_ONLY"] for r in manifest if r["study"]==study}
    c=index(read(src/f"{study}_CPM_projection_QC_LOCAL_ONLY.tsv"),"sample_id")
    q=index(read(prep/study/"library_expression_QC_LOCAL_ONLY.tsv"),"library_id_LOCAL_ONLY")
    assert set(c)==set(q)==ids and len(ids)==n
    if study=="SRP157974":assert set(F)==ids
    for alt in ("TMM_CPM","UQ_CPM"):
        a=index(read(src/f"{study}_{alt}_projection_QC_LOCAL_ONLY.tsv"),"sample_id")
        s=index(read(src/f"{study}_CPM_vs_{alt}_LOCAL_ONLY.tsv"),"sample_id")
        assert set(a)==set(s)==ids
        results=[]
        for k in sorted(ids):
            r=s[k];assert truth(r["comparison_evaluable"])
            f=[num(c[k]["fraction_"+cmp]) for cmp in ["NonSE",*["SE"+str(i) for i in range(1,10)]]]
            g=[num(a[k]["fraction_"+cmp]) for cmp in ["NonSE",*["SE"+str(i) for i in range(1,10)]]]
            assert abs(sum(f)-1)<1e-8 and abs(sum(g)-1)<1e-8
            tvd=num(r["total_variation_distance"])
            assert abs(tvd-sum(abs(x-y) for x,y in zip(f,g))/2)<1e-8
            margin=num(c[k]["top1_minus_top2"])
            assert abs(margin-(sorted(f,reverse=True)[0]-sorted(f,reverse=True)[1]))<1e-8
            same=truth(r["dominant_component_concordant"])
            assert same==(c[k]["dominant_component"]==a[k]["dominant_component"])
            results.append({"tvd":tvd,"margin":margin,"same":same,"k":k,
                            "kl":num(c[k]["normalized_KL_per_input_mass"]),
                            "top10":num(q[k]["raw_top10_fraction"])})
        print("STUDY",study,alt,desc(results))
        if study=="SRP157974":
            groups=defaultdict(list)
            for r in results:groups[F[r["k"]]["author_sequencing_batch_at_base_level"]].append(r)
            assert {k:len(v) for k,v in groups.items()}=={"batch1":180,"batch2":42,"batch3":112}
            for batch,rr in sorted(groups.items()):
                print("FUSCC_BATCH",alt,batch,desc(rr),"CPM_KL",stats([r["kl"] for r in rr]),"TOP10",stats([r["top10"] for r in rr]))
            ok=[r for r in results if F[r["k"]]["rep_suffix_status"]=="AGREES"]
            assert len(ok)==331
            print("FUSCC_SUFFIX_AGREEMENT_SENSITIVITY",alt,desc(ok))
print("V8_EXISTING_PROJECTION_ROBUSTNESS_REVIEW_COMPLETE")
print("DESCRIPTIVE_ONLY_NO_FITS_NO_EXCLUSIONS_NO_CLINICAL_ADMISSIONS")
