#!/usr/bin/env python3
"""Build additive private grouping evidence and targeted unresolved review queue."""
import csv, gzip, hashlib, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert hashlib.sha256(ledger.read_bytes()).hexdigest()=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
with ledger.open(encoding="utf-8",newline="") as f:
    rows=list(csv.DictReader(f,delimiter="\t"))
assert len(rows)==len({r["library_id_LOCAL_ONLY"] for r in rows})==2055
assert all(r["final_approval"]=="NO_UNSIGNED" for r in rows)
remaining=[r for r in rows if r["proposed_non_tumour_exclusion"]=="NO"]
assert len(remaining)==2033
groups=defaultdict(list)
for r in remaining:
    groups[(r["study"],r["historical_patient_key_LOCAL_ONLY"])].append(r)
reused={k:v for k,v in groups.items() if len(v)>1}
assert len(groups)==1965 and len(reused)==36 and sum(map(len,reused.values()))==104
studies=sorted({s for s,_ in reused})
assert len(studies)==11
def norm(x):return str(x or "").strip().casefold()
def present(x):return norm(x) not in {"","na","n/a","none","null","nan","unknown","not available","not applicable","--"}
def counts(xs):return dict(sorted(Counter(xs).items()))
def safe_code(x):
    x=str(x or "").strip()
    return x if re.fullmatch(r"[A-Z][A-Z0-9_]{0,120}",x) else "EMPTY" if not x else "OTHER_SUPPRESSED"
def compare_partition(rs, labels):
    # Compare exact partition of all covered records, not literal key strings.
    forward=defaultdict(set);reverse=defaultdict(set)
    for r in rs:
        e=r["library_id_LOCAL_ONLY"]
        key=r["historical_patient_key_LOCAL_ONLY"]
        if e in labels:
            forward[key].add(labels[e]);reverse[labels[e]].add(key)
    return {"covered":len(labels),"provisional_keys_split_by_source":sum(len(v)>1 for v in forward.values()),
            "source_keys_merging_provisional_keys":sum(len(v)>1 for v in reverse.values()),
            "complete_partition_agreement":len(labels)==len(rs) and all(len(v)==1 for v in forward.values()) and all(len(v)==1 for v in reverse.values())}


# Reviewed scope from the executed source audit; do not infer from labels alone.
expected={
    "SRP139591":(45,12,45,"8dd1a1093e37bdf979c454f1c1d66dddabdb358b41bc767eebdc421c57ff6289"),
    "SRP152822":(21,7,21,"1430f0d78ec66d3ce172dcbe42ab3160906bdf1f6f0f5e4d708fc56f6922c38f"),
    "SRP119676":(29,3,9,"b8ff5f58f81a3a0e5ef523e50b25632246b8fa1d43bc004c59cb07b7f9bd0298"),
    "SRP058107":(8,4,8,"550701c2aaa4937ed06f58268a5ced3b66106db9fe1a73fd4d61c127f0f11c4c"),
    "SRP089805":(7,3,6,"b63664abde485f7358178c99858211755b14cfdc4901e9ecb2708ec830d6ccfc"),
    "SRP036848":(83,2,4,"f7ed0febdf129113ff9e840e684c8333f053cc78f189f36455f66a3c5c2a31d7"),
    "SRP199749":(3,1,3,"146d90ea6ee0b0de779c896c0d641f75e62794913fc5f3a994883ddb36d18dae"),
    "SRP045225":(54,1,2,"e28a97c2aa67a6ca3a1fa2f41627d84319ddb73eb7ac6779eb8d964d0a728dbd"),
    "SRP131744":(3,1,2,"6a646c2983f435ed76377e64eeac3c3f22a3fdb0d74a0874d256dab2e8d7aa0e"),
    "SRP134188":(8,1,2,"1faa9894385d52b1829170ec109b9ee495a64b4b1dabe8be6feda755346c02ec"),
    "SRP163173":(85,1,2,"a1ff0407ae7b252cab577b560e23d4e41efbe3f0ca02ba5bf41753a4cf972986"),
}
supported={
    "SRP139591":("DIRECT_GEO_PATIENT_FIELD_EXACT_SRX","SOURCE_SUPPORTED_EXPLICIT_PATIENT_PARTITION","GSE112996","DISTINCT_REGION_LABELS_PRESERVE_PATIENT_CLUSTER"),
    "SRP152822":("DIRECT_GEO_SAMPLE_TITLE_AND_COHORT_PATIENT_DESIGN","SOURCE_SUPPORTED_PATIENT_REGION_DESIGN","GSE116782","DISTINCT_REGION_LABELS_PRESERVE_PATIENT_CLUSTER"),
    "SRP119676":("DIRECT_GEO_TUMOUR_TISSUE_AND_PATIENT_PREFIX_REPLICATES","SOURCE_SUPPORTED_PATIENT_PREFIX_PARTITION","GSE104766","REGION_VS_TECHNICAL_REPLICATE_REVIEW_REQUIRED"),
}
# Only source-supported grouping is attached. No clinical approval, no collapse.
assert set(studies)==set(expected)
assert len(sys.argv)==3
out=Path(sys.argv[2]).resolve()
assert out.parent==(repo/"private_data").resolve()
assert out.name=="v8_repeated_group_evidence_addendum_v1" and not out.exists()
def attrs(r):
    d=defaultdict(set)
    for seg in r.get("sample_attributes","").split("|"):
        if ";;" not in seg:continue
        k,v=seg.split(";;",1)
        key=re.sub(r"[^a-z0-9]+","_",norm(k)).strip("_")
        if key and present(v):d[key].add(norm(v))
    return d
addendum=[];queue=[];summary=[]
for s in sorted(expected):
    n,ng,nr,digest=expected[s]
    rs=[r for r in remaining if r["study"]==s]
    rgs={k:v for k,v in reused.items() if k[0]==s}
    rr=[r for v in rgs.values() for r in v]
    assert (len(rs),len(rgs),len(rr))==(n,ng,nr)
    path=rec/"raw_studies"/s/"metadata"/f"sra.sra.{s}.MD.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    wanted={r["library_id_LOCAL_ONLY"] for r in rs}
    matched=defaultdict(list)
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r["experiment_acc"] in wanted:matched[r["experiment_acc"]].append(r)
    assert set(matched)==wanted and all(len(v)==1 for v in matched.values())
    native={k:v[0] for k,v in matched.items()}
    av={k:attrs(r) for k,r in native.items()}
    v5={r["library_id_LOCAL_ONLY"]:json.loads(r["v5_original_JSON_LOCAL_ONLY"]) for r in rs}
    if s in supported:
        grade,status,ref,action=supported[s]
        labels={k:norm(v["v5_patient_group_within_study"]) for k,v in v5.items() if present(v["v5_patient_group_within_study"])}
        assert compare_partition(rs,labels)["complete_partition_agreement"]
        assert all(v5[r["library_id_LOCAL_ONLY"]]["v5_patient_group_evidence_grade"]==grade for r in rr)
        for group in rgs.values():
            reg=[norm(v5[r["library_id_LOCAL_ONLY"]]["v5_specimen_region_or_replicate"]) for r in group]
            assert all(map(present,reg)) and len(set(reg))==len(reg)
        if s=="SRP139591":
            labs={k:next(iter(v["patient"])) for k,v in av.items() if len(v.get("patient",set()))==1}
            assert compare_partition(rs,labs)["complete_partition_agreement"]
        selected=rs
    else:
        selected=rr
        status="HOLD_UNRESOLVED_PROVISIONAL_GROUP"
        action="REVIEW_ORIGINAL_SOURCE_GROUP_AND_SPECIMEN_RELATIONSHIP"
        ref="NATIVE_SRA"
        grade=""
        if s=="SRP163173":
            labs={k:next(iter(v["patient_id"])) for k,v in av.items() if len(v.get("patient_id",set()))==1}
            result=compare_partition(rs,labs)
            assert result["covered"]==85 and result["provisional_keys_split_by_source"]==1 and result["source_keys_merging_provisional_keys"]==0
            assert len({labs[r["library_id_LOCAL_ONLY"]] for r in rr})==2
            status="HOLD_SOURCE_PATIENT_ID_VERSUS_PROVISIONAL_GROUP_CONFLICT"
            action="DO_NOT_SPLIT_OR_MERGE_PENDING_AUTHORITATIVE_IDENTITY_EVIDENCE"
            ref="GSE120741_AND_NATIVE_SRA"
    for r in selected:
        k=r["library_id_LOCAL_ONLY"]
        record={"study":s,"library_id_LOCAL_ONLY":k,
                "historical_patient_key_LOCAL_ONLY":r["historical_patient_key_LOCAL_ONLY"],
                "remaining_group_size":r["remaining_provisional_group_size"],
                "group_evidence_status":status,"source_reference":ref,"native_source_sha256":digest,
                "v5_group_evidence_grade":v5[k].get("v5_patient_group_evidence_grade",""),
                "v5_source_group_LOCAL_ONLY":v5[k].get("v5_patient_group_within_study",""),
                "v5_region_or_replicate_LOCAL_ONLY":v5[k].get("v5_specimen_region_or_replicate",""),
                "proposed_group_handling_NOT_APPROVAL":action,
                "patient_independence_certified":"NO","final_clinical_approval":"NO_UNSIGNED"}
        addendum.append(record)
        if s not in supported:
            # Private reviewer queue, no individual values in stdout.
            queue.append({**record,"native_source_record_JSON_LOCAL_ONLY":json.dumps(native[k],ensure_ascii=False,sort_keys=True),
                          "historical_v5_record_JSON_LOCAL_ONLY":json.dumps(v5[k],ensure_ascii=False,sort_keys=True),
                          "reviewer_group_decision":"UNSIGNED","reviewer_evidence_pointer_LOCAL_ONLY":""})
    # Show how well plausible source identifier fields track the full study partition.
    # None of these equality comparisons alone establishes patient identity.
    if s not in supported:
        print("STUDY",s,"UNRESOLVED_REPEATED_GROUPS",ng,"RECORDS",nr)
        for field in ("isolate","cgsr_id","external_id","tag","patient_id"):
            labs={k:next(iter(v[field])) for k,v in av.items() if len(v.get(field,set()))==1}
            if labs:print("CANDIDATE_IDENTIFIER_NOT_CERTIFICATION",field,compare_partition(rs,labs))
        for field in ("library_selection","library_layout","platform_model","library_construction_protocol"):
            vals=[len({norm(native[r["library_id_LOCAL_ONLY"]].get(field,"")) for r in g}) for g in rgs.values()]
            print("WITHIN_GROUP_DISTINCT_VALUES",field,counts(vals))
    summary.append({"study":s,"study_remaining_n":n,"repeated_groups":ng,"repeated_records":nr,
                    "source_supported_group_records":len(selected) if s in supported else 0,
                    "unresolved_repeated_records":0 if s in supported else nr,
                    "source_conflict_records":nr if s=="SRP163173" else 0,"clinical_admissions":0})
assert len(addendum)==124 and len(queue)==29
assert len({r["library_id_LOCAL_ONLY"] for r in addendum})==124
assert sum(r["source_supported_group_records"] for r in summary)==95
assert sum(r["unresolved_repeated_records"] for r in summary)==29
assert sum(r["source_conflict_records"] for r in summary)==2
out.mkdir(mode=0o700)
def write(name,records):
    path=out/name
    with path.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(records)
    path.chmod(0o600)
    return hashlib.sha256(path.read_bytes()).hexdigest()
print("ADDENDUM_SHA256",write("124_group_evidence_addendum_LOCAL_ONLY.tsv",addendum))
print("REVIEW_QUEUE_SHA256",write("29_unresolved_group_review_LOCAL_ONLY.tsv",queue))
print("AGGREGATE_SHA256",write("11_study_group_review_UNSIGNED.tsv",summary))
print("\t".join(summary[0]))
for r in summary:print("\t".join(str(v) for v in r.values()))
print("V8_REPEATED_GROUP_EVIDENCE_ADDENDUM_COMPLETE")
print("NO_FROZEN_INPUT_CHANGES_NO_GROUP_COLLAPSE_NO_CLINICAL_ADMISSIONS")
