#!/usr/bin/env python3
"""One cohort-wide exact accession screen of the frozen administrative 2055; not patient certification."""

import csv
import io
import zipfile
import gzip
import hashlib
import re
import sys
from collections import defaultdict, Counter
from pathlib import Path
import pyarrow.parquet as pq

root = Path(sys.argv[1])
repo = root / 'RECOUNT3_CERTIFICATION_GIT'

files = {
 'atlas': root/'14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz',
 's04': root/'04_candidate_evidence_v1_0_4/04_candidate_sample_full_metadata.parquet',
 's05': root/'05_study_sample_adjudication_v1_0_2_1/18_complete_adjudication_manifest.parquet',
 'cand': root/'41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/03_candidate_integrated_evidence_3602.tsv',
 'full': root/'41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/01_full_source_experiment_identity_5221.tsv',
}

expected = {
 'atlas': '7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092',
 's04': '2e74aabd85a62c3f80835691a2ffa0c11f374bf262142349f579fd4e6f4fc36a',
 's05': '736e89b2969693f504d298c8e7274a1bcf5773a0f18667b114953c575ffd4a32',
 'cand': 'aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81',
 'full': '95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa',
}

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

for name, path in files.items():
    assert path.is_file() and digest(path) == expected[name], (
        'FROZEN_INPUT_MISMATCH_' + name
    )

print('V8_PRIMARY2055_ATLAS_EXACT_SCREEN')
print('FROZEN_INPUT_HASHES_VERIFIED', len(files))

EXP = re.compile(
    r'(?<![A-Za-z0-9])(?:SRX|ERX|DRX)\d+(?![A-Za-z0-9])', re.I
)
BIO = re.compile(
    r'(?<![A-Za-z0-9])(?:SAMN|SAMEA|SAMD)\d+(?![A-Za-z0-9])', re.I
)

def tokens(pattern, value):
    return set(
        x.upper() for x in pattern.findall(str(value or ''))
    )

def read_tsv(path, zipped=False):
    opener = gzip.open if zipped else open
    with opener(str(path), 'rt', encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

atlas_rows = read_tsv(files['atlas'], True)
atlas = {x['sample_id'].strip().upper() for x in atlas_rows}
assert len(atlas_rows) == len(atlas) == 9465


cand=read_tsv(files["cand"])
full=read_tsv(files["full"])
assert (len(cand),len(full))==(3602,5221)
ci={(r["study"],r["library_id"]):r for r in cand}
fi={(r["study"],r["experiment_acc"]):r for r in full}
assert len(ci)==len(cand) and len(fi)==len(full)
v5=root/"26_recount3_external_primary_v1/RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
assert v5.is_file()
print("V5_ZIP_SHA256",digest(v5))
with zipfile.ZipFile(v5) as z:
    names=[n for n in z.namelist() if n.endswith("/results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv")]
    assert len(names)==1
    rows=list(csv.DictReader(io.StringIO(z.read(names[0]).decode("utf-8-sig")),delimiter="\t"))
primary={(r["study"],r["library_id"]):r for r in rows}
assert len(rows)==len(primary)==2055
assert len({s for s,e in primary})==77
assert primary.keys()<=ci.keys() and primary.keys()<=fi.keys()
assert all(ci[k]["technical_hold"]=="NO" for k in primary)
sarcoma={e for s,e in primary if s=="SRP057793"}
assert len(sarcoma)==140
print("PRIMARY_SCREEN_EXACT_JOIN",len(primary),"STUDIES",77)
def graph_from_parquet(path, key):
    pf = pq.ParquetFile(str(path))
    assert key in pf.schema_arrow.names
    assert 'relation' in pf.schema_arrow.names

    seen = {}
    n = duplicate = conflict = 0

    for batch in pf.iter_batches(
        batch_size=8192, columns=[key, 'relation']
    ):
        data = batch.to_pydict()

        for acc, rel in zip(data[key], data['relation']):
            if isinstance(acc, bytes):
                acc = acc.decode('utf-8', 'replace')

            gsm = str(acc or '').strip().upper()

            if gsm not in atlas:
                continue

            n += 1
            pair = (
                frozenset(tokens(EXP, rel)),
                frozenset(tokens(BIO, rel))
            )

            if gsm in seen:
                duplicate += 1
                if pair != seen[gsm]:
                    conflict += 1
            else:
                seen[gsm] = pair

    return seen, n, duplicate, conflict

graphs = []

for label, field in [
    ('s04', 'geo_accession'),
    ('s05', 'sample_id')
]:
    g, n, d, c = graph_from_parquet(files[label], field)

    print(
        'ATLAS_SOURCE', label,
        'MATCHING_ROWS', n,
        'DISTINCT_GSMS', len(g),
        'DUPLICATE_GSM_ROWS', d,
        'CONFLICTS_WITHIN_SOURCE', c,
        'WITH_EXPERIMENT', sum(bool(x[0]) for x in g.values()),
        'WITH_BIOSAMPLE', sum(bool(x[1]) for x in g.values()),
        'WITH_BOTH', sum(bool(x[0] and x[1]) for x in g.values())
    )
    graphs.append((g, n, d, c))

g0, n0, d0, c0 = graphs[0]
g1, n1, d1, c1 = graphs[1]

disagreements = sum(
    g0[gsm] != g1[gsm]
    for gsm in atlas if gsm in g0 and gsm in g1
)

print('PARSED_STAGE04_STAGE05_DISAGREEMENTS', disagreements)

complete = (
    len(g0) == len(g1) == 9465
    and n0 == n1 == 9465
    and d0 == d1 == c0 == c1 == disagreements == 0
    and all(x[0] and x[1] for x in g0.values())
)

print('FULL_EXPERIMENT_AND_BIOSAMPLE_GRAPH_COVERAGE', complete)


assert complete, "INCOMPLETE_ATLAS_GRAPH"
exp_index=defaultdict(set);bio_index=defaultdict(set)
for gsm,(exps,bios) in g0.items():
    for x in exps:exp_index[x].add(gsm)
    for x in bios:bio_index[x].add(gsm)
print("ATLAS_UNIQUE_EXPERIMENT_IDS",len(exp_index))
print("ATLAS_UNIQUE_BIOSAMPLE_IDS",len(bio_index))
def parse(path,series):
    samples={}; cur=None; seen_series=set()
    with gzip.open(path,"rt",encoding="utf-8",errors="strict") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith("^"):
                cur=None
                if line.startswith("^SERIES = "): seen_series.add(line.split(" = ",1)[1])
                if line.startswith("^SAMPLE = "):
                    gsm=line.split(" = ",1)[1]
                    assert re.fullmatch(r"GSM\d+",gsm) and gsm not in samples
                    cur=defaultdict(list);samples[gsm]=cur
            elif cur is not None and line.startswith("!Sample_") and " = " in line:
                k,v=line.split(" = ",1);cur[k].append(v.strip())
    assert series in seen_series and samples, "Incorrect or empty GEO family source"
    return samples
def relations(r):
    return set(re.findall(r"\bSRX\d+\b"," ".join(r.get("!Sample_relation",[]))))

# Two already-locked sarcoma GEO RNA-seq files: add exact GSM alias screen.
cache=repo/"private_data/v8_srp057793_geo_source_v1"
expected_geo={"GSE71119":"249d5e0a299bf35a93802ca7e3a61e2aa8a36f100930300ac9c87400eae87858","GSE75885":"1c59dee75da58ddcff2bc2761f395e3ed27dd5fde83055cc11835dadcb04e6f4"}
sarcoma_gsm=defaultdict(set)
for series,expected_hash in expected_geo.items():
    path=cache/f"{series}_family.soft.gz"
    assert digest(path)==expected_hash
    for gsm,r in parse(path,series).items():
        sr=relations(r)
        if len(sr)==1 and next(iter(sr)) in sarcoma:
            sarcoma_gsm[next(iter(sr))].add(gsm)
assert len(sarcoma_gsm)==138
print("SARCOMA_EXACT_GEO_ALIAS_COVERAGE",len(sarcoma_gsm))
print("SARCOMA_GEO_GSM_ATLAS_HITS",sum(bool(gs&atlas) for gs in sarcoma_gsm.values()))
stats=defaultdict(Counter)
total=Counter()
for k,r in primary.items():
    study,exp=k
    f=fi[k];c=ci[k]
    bfull=tokens(BIO,f["biosamples_all_sources"])
    bcand=tokens(BIO,c["biosample_IDs"])
    bv5=tokens(BIO,r["biosample"])
    ef=bool(exp_index.get(exp.upper()))
    bf=any(bio_index.get(x) for x in bfull)
    bc=any(bio_index.get(x) for x in bcand)
    bv=any(bio_index.get(x) for x in bv5)
    gh=study=="SRP057793" and bool(sarcoma_gsm.get(exp,set())&atlas)
    flags={"N":True,"EXPERIMENT_HIT":ef,"FULL_BIO_HIT":bf,"CAND_BIO_HIT":bc,"V5_BIO_HIT":bv,"SARCOMA_GEO_GSM_HIT":gh,"ANY_EXACT_HIT":ef or bf or bc or bv or gh,"FULL_CAND_BIO_DISCORD":bfull!=bcand,"FULL_V5_BIO_DISCORD":bfull!=bv5,"FULL_BIO_EMPTY":not bfull,"CAND_BIO_EMPTY":not bcand,"V5_BIO_EMPTY":not bv5}
    for name,val in flags.items():
        stats[study][name]+=int(val);total[name]+=int(val)
fields=list(total)
print("STUDY_LEVEL_EXACT_SCREEN")
print("\t".join(["study"]+fields))
for s,ct in sorted(stats.items()):
    print("\t".join([s]+[str(ct[f]) for f in fields]))
print("TOTAL",dict(total))
print("STUDIES_WITH_ANY_EXACT_HIT",sum(c["ANY_EXACT_HIT"]>0 for c in stats.values()))
print("STUDIES_WITH_BIOSAMPLE_FIELD_DISCORDANCE",sum(c["FULL_CAND_BIO_DISCORD"]>0 or c["FULL_V5_BIO_DISCORD"]>0 for c in stats.values()))
print("V8_PRIMARY2055_ATLAS_EXACT_SCREEN_COMPLETE")
print("EXACT_ACCESSION_NONOVERLAP_IS_NOT_PATIENT_INDEPENDENCE")
print("NO_FILES_MODIFIED_NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0")
