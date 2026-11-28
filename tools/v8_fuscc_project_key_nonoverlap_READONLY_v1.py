#!/usr/bin/env python3
"""New patient-project token screen; not a proof of global patient independence."""
import csv,gzip,hashlib,re,sys
from collections import Counter,defaultdict
from pathlib import Path
import pyarrow.parquet as pq
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def locked(p,h):
    assert sha(p)==h,"INPUT_HASH_MISMATCH "+str(p)
    return p
cross=locked(repo/"private_data/v8_fuscc_patient_crosswalk_addendum_v1/360_fuscc_patient_clinical_crosswalk_LOCAL_ONLY.tsv","8d34f94ba5834e6dc477fb740547ab46801cc89185e1b724ec49356d50400cd5")
rows=read(cross)
keys={r["clinical_project_id_LOCAL_ONLY"] for r in rows}
assert len(rows)==len(keys)==360
assert all(re.fullmatch(r"FUSCCTNBC\d+",k) for k in keys)
assert all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in rows)
disposition={r["clinical_project_id_LOCAL_ONLY"]:r["proposed_review_disposition"] for r in rows}
assert Counter(disposition.values())=={"PROPOSED_PATIENT_LINKED_REVIEW_CANDIDATE":331,"HOLD_LIBRARY_SUFFIX_PROVENANCE":3,"HISTORICAL_QC_HOLD_PRESERVED":26}
atlasfile=locked(root/"14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz","7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092")
atlasrows=read(atlasfile)
atlas={r["sample_id"].strip().upper() for r in atlasrows}
assert len(atlasrows)==len(atlas)==9465
PAT=re.compile(r"(?<![A-Za-z0-9])FUSCCTNBC\d+(?![A-Za-z0-9])",re.I)
def tokens(x):
    # Nested values rendered as text; only exact namespaced project tokens.
    return {t.upper() for t in PAT.findall(str(x or ""))}
def scan(rows,label,idfield,restrict=None):
    hits=defaultdict(set);other=defaultdict(set);fields=Counter();seen=Counter();n=0
    for r in rows:
        identifier=str(r[idfield] or "").strip().upper()
        if restrict is not None and identifier not in restrict:continue
        seen[identifier]+=1;n+=1
        for f,v in r.items():
            ts=tokens(v)
            if ts&keys:
                fields[f]+=1
                for k in ts&keys:hits[k].add(identifier)
            for k in ts-keys:other[k].add(identifier)
    print("SOURCE",label,"ROWS",n,"DISTINCT_IDS",len(seen),"DUPLICATE_ROWS",n-len(seen),"EXACT_PROJECT_KEYS_HIT",len(hits),"HIT_ROWS",len(set().union(*hits.values())) if hits else 0,"OTHER_FUSCCTNBC_KEYS",len(other))
    print("MATCHED_FIELDS_ROW_COUNTS",label,dict(fields))
    print("HIT_DISPOSITIONS",label,dict(Counter(disposition[k] for k in hits)))
    return seen,hits
scan(atlasrows,"ATLAS_CORE","sample_id",atlas)
sources=[
("STAGE04",root/"04_candidate_evidence_v1_0_4/04_candidate_sample_full_metadata.parquet","2e74aabd85a62c3f80835691a2ffa0c11f374bf262142349f579fd4e6f4fc36a","geo_accession"),
("STAGE05",root/"05_study_sample_adjudication_v1_0_2_1/18_complete_adjudication_manifest.parquet","736e89b2969693f504d298c8e7274a1bcf5773a0f18667b114953c575ffd4a32","sample_id")]
allhits=set()
for label,path,h,idfield in sources:
    pf=pq.ParquetFile(locked(path,h))
    assert idfield in pf.schema_arrow.names
    print("METADATA_COLUMNS_SCANNED",label,len(pf.schema_arrow.names))
    def iterrows():
        for batch in pf.iter_batches(batch_size=1024):
            yield from batch.to_pylist()
    seen,hits=scan(iterrows(),label,idfield,atlas)
    assert set(seen)==atlas and all(n==1 for n in seen.values()),"FROZEN_METADATA_AXIS"
    allhits.update(hits)
ledger=locked(repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv","943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0")
lr=read(ledger);assert len(lr)==2055
others=[r for r in lr if r["study"]!="SRP157974"]
assert len(others)==1721
seen,hits=scan(others,"OTHER_EXTERNAL_PRIMARY_SCREEN","library_id_LOCAL_ONLY")
assert len(seen)==1721
print("ATLAS_METADATA_UNION_PROJECT_KEYS_HIT",len(allhits))
print("OTHER_EXTERNAL_PROJECT_KEYS_HIT",len(hits))
print("V8_FUSCC_PROJECT_KEY_NONOVERLAP_REVIEW_COMPLETE")
print("EXACT_TOKEN_SCREEN_NOT_GLOBAL_PATIENT_INDEPENDENCE_NO_CLINICAL_ADMISSIONS")
