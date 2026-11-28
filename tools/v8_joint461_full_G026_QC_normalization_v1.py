#!/usr/bin/env python3
"""Four-study joint461 technical preparation, NOT a clinical or projection release."""
import csv,gzip,json,hashlib,re,sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from v7r_prepare_configured_matched_G026_full499_v1 import (
    digest,read_tsv,original_module,record_matrix,EXPECTED_MAP_SHA,EXPECTED_STAGE09_SHA)
STUDIES={"SRP157974":(334,334,448),"SRP057793":(90,90,149),"SRP002628":(14,14,30),"SRP094795":(23,45,47)}
RUN=re.compile(r"(?<![A-Za-z0-9])(?:SRR|ERR|DRR)\d+(?![A-Za-z0-9])",re.I)
READY="INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED"
HOLD="HOLD_INPUT_CONTRACT_UNRESOLVED"
def require(ok,message):
    if not ok:raise ValueError(message)
def check(p,h):
    require(bool(re.fullmatch(r"[0-9a-f]{64}",h)),"INVALID_EXPECTED_SHA "+str(p))
    observed=digest(p)
    require(observed==h,"HASH_MISMATCH "+str(p)+" expected="+h+" observed="+observed)
def from_log(log,prefix):
    hits=re.findall(r"^"+re.escape(prefix)+r" ([0-9a-f]{64})$",log,re.M)
    require(len(hits)==1,"MISSING_OR_REPEATED_LOG_DIGEST "+prefix)
    return hits[0]
def verified_audit(path,marker,required):
    # Only trust the existing checksum manifest after checking the actual files.
    entries={}
    for line in (path/"SHA256SUMS.txt").read_text().splitlines():
        m=re.fullmatch(r"([0-9a-f]{64}) [ *](.+)",line)
        require(m is not None,"BAD_AUDIT_CHECKSUM_LINE")
        fp=Path(m[2]);fp=(path/fp).resolve() if not fp.is_absolute() else fp.resolve()
        require(fp not in entries,"DUPLICATE_AUDIT_PATH")
        entries[fp]=m[1]
    require({x.resolve() for x in required}<=entries.keys(),"REQUIRED_AUDIT_PATH_MISSING "+str(path))
    for fp,h in entries.items():check(fp,h)
    log=(path/"stdout.txt").read_text()
    require(marker in log.splitlines(),"AUDIT_MARKER_MISSING "+marker)
    return log
def write_tsv(path,rows):
    require(bool(rows),"EMPTY_OUTPUT "+str(path))
    with path.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    path.chmod(0o600)
def stats(x):
    a=np.asarray(x,dtype=float)
    require(a.size>0 and np.isfinite(a).all(),"INVALID_QC_VECTOR")
    return {"min":float(a.min()),"median":float(np.median(a)),"max":float(a.max())}
def load_native(p):
    with gzip.open(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def run_tokens(v):return {x.upper() for x in RUN.findall(v or "")}
def amend(rows,native):
    grouped=defaultdict(list)
    for r in native:grouped[r["experiment_acc"]].append(r)
    fields=("study","sample_acc","experiment_acc","submission_acc","submission_center",
            "experiment_title","library_strategy","library_source","library_selection","library_layout",
            "library_construction_protocol","platform_model","sample_attributes","experiment_attributes","sample_name","sample_title")
    amended=[];n=0
    for src in rows:
        r=dict(src);r["original_input_status"]=r["input_status"];r["cardinality_evidence"]="ORIGINAL_SINGLE_NATIVE_ROW"
        if r["study"]=="SRP094795":
            a=grouped[r["library_id_LOCAL_ONLY"]]
            expected=set(r["stage2_run_ids_LOCAL_ONLY"].split("|"))
            per=[run_tokens(x["external_id"]) for x in a]
            require(len(a)==int(r["stage2_expected_runs"])==int(r["stage2_verified_runs"])==int(r["native_exact_experiment_rows"])==len(expected),"CERVICAL_NATIVE_COUNT_CONFLICT")
            require(all(len(x)==1 for x in per) and set().union(*per)==expected,"CERVICAL_NATIVE_RUN_SET_CONFLICT")
            require(all(len({x[f].strip() for x in a})==1 and a[0][f].strip() for f in fields),"CERVICAL_SPECIMEN_OR_LIBRARY_FIELD_CONFLICT")
            require(all(x["study"]=="SRP094795" and x["experiment_acc"]==r["library_id_LOCAL_ONLY"] for x in a),"CERVICAL_EXPERIMENT_CONFLICT")
            require(r["run_membership_and_count_ok"]==r["raw_expected_runs_present"]=="TRUE" and r["historical_hold_queue_rows"]=="0","OTHER_INPUT_HOLD_PRESENT")
            if r["input_status"]==HOLD:
                require(len(a)==2,"UNEXPECTED_HOLD_TYPE")
                r["input_status"]=READY
                r["native_library_selection"]=a[0]["library_selection"]
                r["native_platform"]=a[0]["platform_model"]
                r["cardinality_evidence"]="EXACT_TWO_RUNS_ONE_EXPERIMENT_INVARIANT_SOURCE"
                n+=1
            else:require(len(a)==1 and r["input_status"]==READY,"UNEXPECTED_SINGLE_RUN_STATUS")
        require(r["input_status"]==READY and r["final_clinical_approval"]=="NO_UNSIGNED" and r["patient_independence_certified"]=="NO","UNRESOLVED_OR_CHANGED_CLINICAL_STATUS")
        amended.append(r)
    require(n==22,"UNEXPECTED_AMENDMENT_COUNT")
    return amended
def main():
    require(len(sys.argv)==3,"ROOT OUTPUT required")
    root=Path(sys.argv[1]).resolve();repo=root/"RECOUNT3_CERTIFICATION_GIT";rec=root/"26_recount3_external_primary_v1"
    out=Path(sys.argv[2]).resolve()
    require(out.parent==(repo/"private_data").resolve() and out.name=="v8_joint461_full_G026_QC_normalization_v1" and not out.exists(),"OUTPUT_PROTECTION")
    a=repo/"diagnostics/20261128_v8_joint461_technical_input_manifest_v1_0_1"
    mp=repo/"private_data/v8_joint461_technical_input_manifest_v1_0_1/461_technical_input_manifest_LOCAL_ONLY.tsv"
    log=verified_audit(a,"V8_JOINT461_TECHNICAL_INPUT_MANIFEST_V1_0_1_COMPLETE",[mp,a/"stdout.txt",repo/"tools/v8_joint461_technical_input_manifest_v1_0_1.py"])
    ca=repo/"diagnostics/20261128_v8_srp094795_native_run_cardinality_v1"
    verified_audit(ca,"V8_SRP094795_NATIVE_RUN_CARDINALITY_COMPLETE",[ca/"stdout.txt",repo/"tools/v8_srp094795_native_run_cardinality_READONLY_v1.py"])
    _,rows=read_tsv(mp)
    require(len(rows)==len({r["library_id_LOCAL_ONLY"] for r in rows})==461,"MANIFEST_IDENTITY")
    require(Counter(r["study"] for r in rows)=={s:x[0] for s,x in STUDIES.items()},"STUDY_COUNTS")
    require(Counter(r["input_status"] for r in rows)=={READY:439,HOLD:22},"ORIGINAL_STATUS_COUNTS")
    native=rec/"raw_studies/SRP094795/metadata/sra.sra.SRP094795.MD.gz"
    check(native,from_log(log,"NATIVE_SOURCE SRP094795 SHA256"))
    rows=amend(rows,load_native(native))
    require(len({r["library_id_LOCAL_ONLY"] for r in rows})==461,"AMENDED_IDENTITY")
    allruns=[x for r in rows for x in r["stage2_run_ids_LOCAL_ONLY"].split("|")]
    require(len(allruns)==len(set(allruns))==483 and all(RUN.fullmatch(x) for x in allruns),"RUN_IDENTITY_GLOBAL")
    qpath=rec/"stage2_output/run_qc_all.tsv"
    check(qpath,from_log(log,"INPUT_SHA256_PASS run_qc_all.tsv"))
    _,qc=read_tsv(qpath);qmap=defaultdict(list)
    for r in qc:qmap[r["run_id"].strip().upper()].append(r)
    gmap=rec/"stage5_output/00_gene_row_mapping.tsv";check(gmap,EXPECTED_MAP_SHA)
    _,mapping=read_tsv(gmap)
    require(len(mapping)==63856 and [int(r["row_zero_based"]) for r in mapping]==list(range(63856)),"G026_ROW_MAP")
    symbols=[r["gene_symbol"] for r in mapping]
    rel=root/"RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
    frozen=rel/"09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py";check(frozen,EXPECTED_STAGE09_SHA)
    # Lock the gene ledger against the already executed prostate matched preparation.
    v7r=repo/"private_data/v7r_SRP157215_prepared_v1/input_provenance.json"
    v7audit=repo/"diagnostics/20261127_v7r_SRP157215_calibration_v1"
    v7entries=[l.split(maxsplit=1) for l in (v7audit/"SHA256SUMS.txt").read_text().splitlines()]
    matches=[h for h,name in v7entries if Path(name.lstrip("* ")).resolve()==v7r.resolve()]
    require(len(matches)==1,"V7R_PROVENANCE_CHECKSUM_NOT_FOUND")
    check(v7r,matches[0])
    v7=json.loads(v7r.read_text())
    require(v7["study"]=="SRP157215" and v7["n_model_genes"]==499,"V7R_REFERENCE_CONTRACT")
    candidates=[root/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv",rel/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv"]
    geneaudit=next((x for x in candidates if x.is_file()),None);require(geneaudit is not None,"MODEL_GENE_LEDGER_MISSING")
    check(geneaudit,v7["frozen_model_gene_audit_sha256"])
    _,avail=read_tsv(geneaudit);genes=[r["model_gene"].strip() for r in avail]
    require(len(genes)==len(set(genes))==499,"MODEL_GENE_IDENTITY")
    old=original_module(frozen)
    out.mkdir(mode=0o700)
    write_tsv(out/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv",rows)
    print("CERVICAL_CARDINALITY_AMENDMENT_PASS",22,"TOTAL",461,"RUNS",483,flush=True)
    aggregate=[];provenance={}
    for study,(n,nrun,nraw) in STUDIES.items():
        rs=sorted((r for r in rows if r["study"]==study),key=lambda r:r["library_id_LOCAL_ONLY"])
        ids=[r["library_id_LOCAL_ONLY"] for r in rs];owner={}
        for i,r in enumerate(rs):
            runs=r["stage2_run_ids_LOCAL_ONLY"].split("|")
            require(len(runs)==int(r["stage2_expected_runs"])==int(r["stage2_verified_runs"]),"RUN_COUNT")
            for run in runs:
                require(run not in owner and len(qmap[run])==1,"STAGE2_RUN_DUPLICATE_OR_MISSING")
                q=qmap[run][0]
                require(q["study"]==study and q["library_id"]==ids[i],"STAGE2_RUN_OWNER_CONFLICT")
                owner[run]=i
        require(len(ids)==n and len(owner)==nrun,"STUDY_MEMBERSHIP")
        nativepath=rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz"
        check(nativepath,from_log(log,"NATIVE_SOURCE "+study+" SHA256"))
        raw=rec/"raw_studies"/study/"gene_sums"/f"sra.gene_sums.{study}.G026.gz"
        rawsha=digest(raw)
        rawmat=np.zeros((63856,n),dtype=np.float64)
        sums=np.zeros(nrun,dtype=np.float64);det=np.zeros(nrun,dtype=np.int64)
        gene_order_version_differences=0;max_mass_error=0.
        with gzip.open(raw,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
            header=None
            for _ in range(30):
                line=f.readline()
                if not line:break
                parts=line.rstrip("\r\n").split("\t")
                if len(parts)==nraw+1 and all(RUN.fullmatch(x.strip().strip('"')) for x in parts[1:]):
                    header=[x.strip().strip('"').upper() for x in parts[1:]]
                    headerhash=hashlib.sha256(line.encode()).hexdigest()
                    break
            require(header is not None and len(set(header))==nraw,"RAW_HEADER_CONTRACT "+study)
            require(all(r["raw_header_sha256"]==headerhash for r in rs),"RAW_HEADER_CHANGED "+study)
            require(set(owner)<=set(header),"RAW_RUN_MISSING "+study)
            positions=np.array([j for j,x in enumerate(header) if x in owner],dtype=np.int64)
            selected_runs=[header[j] for j in positions]
            groups=np.array([owner[x] for x in selected_runs],dtype=np.int64)
            require(len(positions)==nrun,"SELECTED_RUN_AXIS")
            rowcount=0
            for line in f:
                require(rowcount<63856,"EXTRA_G026_ROWS "+study)
                vals=line.rstrip("\r\n").split("\t")
                require(len(vals)==nraw+1,"G026_COLUMN_COUNT "+study)
                name=vals[0].strip().strip('"');mr=mapping[rowcount]
                require(name==mr["source_gene_id"] or name.split(".",1)[0]==mr["gene_base"],"G026_GENE_AXIS "+study+" row="+str(rowcount))
                gene_order_version_differences+=int(name!=mr["source_gene_id"])
                x=np.asarray(vals[1:],dtype=np.float64)
                require(np.isfinite(x).all() and (x>=0).all(),"G026_INVALID_VALUES "+study)
                y=x[positions];sums+=y;det+=(y>0)
                agg=np.bincount(groups,weights=y,minlength=n)
                rawmat[rowcount,:]=agg
                max_mass_error=max(max_mass_error,abs(float(y.sum()-agg.sum()))/max(1.,float(y.sum())))
                rowcount+=1
            require(rowcount==63856 and max_mass_error<=1e-12,"G026_ROW_OR_MASS_CONTRACT "+study)
        check(raw,rawsha)
        # Exact detection count and relative coverage comparison with historical Stage2 run QC.
        expected_det=np.array([int(qmap[k][0]["n_genes_coverage_gt0"]) for k in selected_runs])
        expected_mass=np.array([float(qmap[k][0]["gene_body_coverage_bp"]) for k in selected_runs])
        require(np.isfinite(expected_mass).all() and (expected_mass>0).all(),"STAGE2_COVERAGE_INVALID")
        relerr=np.abs(sums-expected_mass)/np.maximum(1.,expected_mass)
        bad_det=int(np.count_nonzero(det!=expected_det));bad_mass=int(np.count_nonzero(relerr>1e-8))
        print("RAW_QC",study,"ROWS",rowcount,"RUNS",nrun,"LIBRARIES",n,"DETECTION_MISMATCH",bad_det,"COVERAGE_MISMATCH",bad_mass,"MAX_REL_COVERAGE_ERROR",float(relerr.max()),"AGGREGATION_MASS_ERROR",max_mass_error,"GENE_VERSION_TEXT_DIFFERENCES",gene_order_version_differences,flush=True)
        require(bad_det==bad_mass==0,"RAW_STAGE2_QC_DISAGREEMENT "+study)
        d=out/study;d.mkdir(mode=0o700)
        # Store aggregated all-gene coverage for future genomewide work; no biological normalization of raw.
        np.savez_compressed(d/"63856_G026_library_coverage_LOCAL_ONLY.npz",coverage=rawmat,library_ids=np.asarray(ids),source_gene_ids=np.asarray([r["source_gene_id"] for r in mapping]))
        factors=record_matrix(rawmat,symbols,genes,old,"recount3",ids,d)
        maps=old.build_gene_maps(pd.Index(symbols),genes)
        model=np.zeros((499,n),dtype=float)
        np.add.at(model,maps["model_row_indices"],rawmat[maps["model_row_positions"],:])
        library_mass=rawmat.sum(axis=0);library_detection=(rawmat>0).sum(axis=0)
        top10=np.partition(rawmat,-10,axis=0)[-10:,:].sum(axis=0)/library_mass
        model_detection=(model>0).sum(axis=0)
        qrows=[]
        for i,r in enumerate(rs):
            qrows.append({"study":study,"library_id_LOCAL_ONLY":ids[i],"n_technical_runs":r["stage2_expected_runs"],"raw_gene_rows_detected":int(library_detection[i]),"raw_gene_detection_fraction":float(library_detection[i]/63856),"raw_coverage_bp":float(library_mass[i]),"raw_top10_fraction":float(top10[i]),"valid_symbol_coverage_bp":float(factors["library_mass"][i]),"model_genes_detected":int(model_detection[i]),"TMM_factor":float(factors["tmm"][i]),"UQ_factor":float(factors["uq"][i]),"UQ_evaluable":str(bool(factors["uq_valid"][i])).upper(),"clinical_approval":"NO_UNSIGNED"})
        write_tsv(d/"library_expression_QC_LOCAL_ONLY.tsv",qrows)
        # Descriptive flags only, not thresholds for excluding biological samples.
        summary={"study":study,"libraries":n,"technical_runs":nrun,"raw_run_columns":nraw,"raw_gene_rows":63856,"model_genes_mapped":499,"model_gene_detection_min":int(model_detection.min()),"model_gene_detection_median":float(np.median(model_detection)),"model_gene_detection_max":int(model_detection.max()),"raw_gene_detection_fraction_median":float(np.median(library_detection/63856)),"raw_coverage_bp_median":float(np.median(library_mass)),"raw_top10_fraction_max":float(top10.max()),"UQ_evaluable":int(factors["uq_valid"].sum()),"clinical_admissions":0}
        aggregate.append(summary)
        print("LIBRARY_QC",study,"DETECTION",stats(library_detection/63856),"COVERAGE",stats(library_mass),"TOP10",stats(top10),"MODEL_DETECTION",stats(model_detection),flush=True)
        print("STUDY_PREPARATION_COMPLETE",study,summary,flush=True)
        for f in d.iterdir():f.chmod(0o600)
        provenance[study]={"raw_sha256":rawsha,"native_sha256":digest(nativepath),"raw_header_sha256":headerhash,"TMM_reference_index_LOCAL_ONLY":factors["tmm_ref_index_local_only"],"outputs_sha256":{f.name:digest(f) for f in sorted(d.iterdir()) if f.is_file()}}
    write_tsv(out/"4_study_expression_QC_aggregate_UNSIGNED.tsv",aggregate)
    result={"version":"v8_joint461_full_G026_QC_normalization_v1","membership":461,"technical_runs":483,"original_manifest_sha256":digest(mp),"amended_manifest_sha256":digest(out/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv"),"original_stage09_sha256":digest(frozen),"gene_map_sha256":digest(gmap),"model_gene_ledger_sha256":digest(geneaudit),"stage2_run_qc_sha256":digest(qpath),"normalization_scope":"within study, CPM primary TMM/UQ sensitivity; no refitted anchors","aggregate_sha256":digest(out/"4_study_expression_QC_aggregate_UNSIGNED.tsv"),"studies":provenance,"clinical_admissions":0}
    target=out/"input_provenance.json";target.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n");target.chmod(0o600)
    print("V8_JOINT461_FULL_G026_PREPARATION_COMPLETE_PENDING_REVIEW",flush=True)
    print("NO_MODEL_EXECUTION_NO_FROZEN_INPUT_MODIFICATIONS_NO_CLINICAL_ADMISSIONS",flush=True)
if __name__=="__main__":main()
