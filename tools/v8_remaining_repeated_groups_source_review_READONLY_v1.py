#!/usr/bin/env python3
"""Review remaining repeated provisional groups against existing V5 and native source evidence.
No group collapse, certification, or individual identifiers in stdout.
"""
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

print("V8_REMAINING_REPEATED_GROUPS_SOURCE_REVIEW")
print("LOCKED_LEDGER_PASS",len(rows),"REMAINING",len(remaining),"REUSED_GROUPS",len(reused),"REUSED_RECORDS",104,"STUDIES",11)
for s in sorted(studies,key=lambda s:-sum(len(v) for (study,_),v in reused.items() if study==s)):
    rs=[r for r in remaining if r["study"]==s]
    repeated=[v for (study,_),v in reused.items() if study==s]
    rr=[r for v in repeated for r in v]
    originals={r["library_id_LOCAL_ONLY"]:json.loads(r["v5_original_JSON_LOCAL_ONLY"]) for r in rs}
    assert all(v["library_id"]==k and v["study"]==s for k,v in originals.items())
    print("\nSTUDY",s,"REMAINING",len(rs),"REUSED_GROUPS",len(repeated),"REUSED_RECORDS",len(rr),"MULTIPLICITIES",counts(len(v) for v in repeated))
    labels={k:norm(v["v5_patient_group_within_study"]) for k,v in originals.items() if present(v["v5_patient_group_within_study"])}
    print("V5_SOURCE_GROUP_PARTITION",compare_partition(rs,labels))
    for field in ("v5_patient_group_evidence_grade","v5_specimen_evidence_class","patient_evidence_grade","treatment_evidence_grade"):
        print("REUSED_RECORD_GRADE",field,counts(safe_code(originals[r["library_id_LOCAL_ONLY"]].get(field)) for r in rr))
    region_counts=Counter()
    for v in repeated:
        rlabels=[norm(originals[r["library_id_LOCAL_ONLY"]]["v5_specimen_region_or_replicate"]) for r in v]
        if not all(map(present,rlabels)):region_counts["REGION_LABEL_INCOMPLETE"]+=1
        elif len(set(rlabels))==len(v):region_counts["ALL_DISTINCT_SOURCE_REGION_OR_REPLICATE_LABELS"]+=1
        else:region_counts["REPEATED_REGION_OR_REPLICATE_LABEL"]+=1
    print("V5_REGION_GROUP_STRUCTURE",dict(region_counts))
    path=rec/"raw_studies"/s/"metadata"/f"sra.sra.{s}.MD.gz"
    if not path.is_file():
        print("NATIVE_SOURCE_MISSING");continue
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        source=list(csv.DictReader(f,delimiter="\t"))
    wanted=set(originals)
    matched=defaultdict(list)
    for r in source:
        if r["experiment_acc"] in wanted:matched[r["experiment_acc"]].append(r)
    print("NATIVE_SHA256",hashlib.sha256(path.read_bytes()).hexdigest())
    print("NATIVE_EXACT_MATCHED",len(matched),"MISSING",len(wanted-matched.keys()),"MULTIPLE_METADATA_ROWS",sum(len(v)>1 for v in matched.values()))
    attrs=defaultdict(lambda:defaultdict(set))
    for k,rsrc in matched.items():
        for r in rsrc:
            for segment in r.get("sample_attributes","").split("|"):
                if ";;" not in segment:continue
                name,value=segment.split(";;",1)
                name=re.sub(r"[^a-z0-9]+","_",norm(name)).strip("_")
                if name and present(value):attrs[name][k].add(norm(value))
    patient_fields={"patient","patient_id","patient_identifier","patient_number","patient_no","subject","subject_id","subject_identifier","subject_number","donor","donor_id","individual","individual_id","participant","participant_id"}
    print("NATIVE_ATTRIBUTE_KEYS",sorted(attrs))
    for name,by in sorted(attrs.items()):
        if name not in patient_fields:continue
        labels={k:next(iter(v)) for k,v in by.items() if len(v)==1}
        print("EXPLICIT_NATIVE_SUBJECT_FIELD",name,"MULTIVALUE",sum(len(v)>1 for v in by.values()),compare_partition(rs,labels))
    # Count categorical native assay labels without sample descriptors.
    for field in ("library_strategy","library_selection","library_layout"):
        print("NATIVE_ASSAY_LABEL",field,counts(safe_code(str(r.get(field,"")).upper().replace("-","_").replace(" ","_")) for vals in matched.values() for r in vals))
print("\nV8_REMAINING_REPEATED_GROUPS_SOURCE_REVIEW_COMPLETE")
print("NO_FILES_MODIFIED_NO_GROUP_COLLAPSE_NO_CLINICAL_ADMISSIONS")
