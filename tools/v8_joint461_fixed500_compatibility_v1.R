#!/usr/bin/env Rscript
# Technical investigation only; no clinical release, no matched assay claims.
suppressPackageStartupMessages({library(data.table);library(SpatialEcoTyper);library(NMF);library(jsonlite)})
argv <- commandArgs(trailingOnly=TRUE)
stopifnot(length(argv)==3L)
ROOT <- normalizePath(argv[[1]],mustWork=TRUE)
SRC <- normalizePath(argv[[2]],mustWork=TRUE)
FINAL <- normalizePath(argv[[3]],mustWork=TRUE)
REPO <- file.path(ROOT,"RECOUNT3_CERTIFICATION_GIT")
stopifnot(!identical(SRC,FINAL),length(list.files(FINAL,all.files=TRUE,no..=TRUE))==0L)
# Evaluate only selected function definitions from the checksum-verified V7R script.
# Never source its paired experiment execution statements.
exprs <- as.list(parse(file=file.path(REPO,"tools/v7r_project_configured_matched_full500_KL_v1.R")))
for(nm in c("REL","S9","S10","BRANCHES","COMPONENTS","CLR_FLOOR","RECONSTRUCTION_EPS",
            "SUM_TO_ONE_TOLERANCE","W_FILE","assert","find_frozen","branch_assets","key",
            "read_499","read_params","scaled_model","archiveZ","harmonize_expression_to_W",
            "score_fractions","clr_matrix","dominant_component","calculate_reconstruction_adequacy",
            "normalization_sensitivity","full500_bundle","fixed500_project_raw_h")){
  hit <- Filter(function(e)is.call(e)&&identical(e[[1L]],as.name("<-"))&&length(e)==3L&&identical(e[[2L]],as.name(nm)),exprs)
  stopifnot(length(hit)==1L)
  eval(hit[[1L]],envir=globalenv())
}
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3","Package version mismatch")
W <- readRDS(W_FILE)
assert(identical(W,readRDS(file.path(system.file("extdata",package="SpatialEcoTyper"),"Bulk_SE_Recovery_W.rds"))),"Installed W mismatch")
assert(is.matrix(W)&&identical(dim(W),c(500L,10L))&&identical(colnames(W),COMPONENTS)&&all(is.finite(W))&&all(W>=0),"W contract")
W_genes <- unique(sub("__.*$","",rownames(W)))
assert(length(W_genes)==499L,"W base gene contract")
gene_audit <- fread(find_frozen(S9,"08_model_gene_availability.tsv"),sep="\t",quote="")
frozen_mapping <- fread(find_frozen(S10,"16A_frozen_W_gene_symbol_harmonization.tsv"),sep="\t",quote="")
assert(nrow(gene_audit)==499L&&nrow(frozen_mapping)==1497L,"Frozen mapping dimensions")
prepare_z <- function(X,branch){
  z <- harmonize_expression_to_W(scaled_model(X,read_params(branch)),W,gene_audit,branch)
  ledger <- copy(frozen_mapping[as.character(frozen_mapping[["branch"]])==branch])
  setorder(ledger,W_gene_index)
  assert(nrow(ledger)==499L&&identical(as.integer(ledger$W_gene_index),seq_len(499L))&&
         identical(as.character(ledger$W_model_gene),W_genes)&&
         identical(as.character(z$mapping$input_model_gene),as.character(ledger$input_model_gene)),"Frozen harmonization mismatch")
  assert(identical(rownames(z$matrix),W_genes),"Z gene order")
  z$matrix
}
# Frozen reference replay, same 22 calibration specimens, not new candidates.
control <- fread(file.path(REPO,"private_data/v7r_SRP157215_prepared_v1/pair_manifest_LOCAL_ONLY.tsv"),sep="\t",quote="")
gsm <- as.character(control$ARCHS4_GSM);N_PAIRS <- 22L
assert(length(gsm)==22L&&uniqueN(gsm)==22L,"Calibration specimen count")
for(branch in BRANCHES){
  X <- read_499(find_frozen(S9,branch_assets[[branch]]$expr),gsm,full_atlas=TRUE)
  Z <- prepare_z(X,branch)
  err <- max(abs(Z-archiveZ(branch)))
  assert(is.finite(err)&&err<=1e-8,paste("Frozen anchor replay failure",branch))
  htab <- fread(find_frozen(S10,branch_assets[[branch]]$h),sep="\t",quote="")
  ftab <- fread(find_frozen(S10,branch_assets[[branch]]$f),sep="\t",quote="")
  assert(nrow(htab)==9465L&&nrow(ftab)==9465L&&uniqueN(htab$sample_id)==9465L&&uniqueN(ftab$sample_id)==9465L,"Archive dimensions")
  ix <- match(gsm,htab$sample_id);jx <- match(gsm,ftab$sample_id)
  assert(!anyNA(ix)&&!anyNA(jx),"Calibration archive IDs")
  ah <- as.matrix(htab[ix,..COMPONENTS]);af <- as.matrix(ftab[jx,..COMPONENTS])
  dimnames(ah) <- dimnames(af) <- list(gsm,COMPONENTS)
  assert(max(abs(score_fractions(ah)-af))<=1e-8,"Archived H/fraction disagreement")
  H <- fixed500_project_raw_h(Z,W,paste0("ARCHS4_CONTROL_",branch),h_init=ah)
  tvd <- max(rowSums(abs(score_fractions(H)-af))/2)
  assert(is.finite(tvd)&&tvd<=5e-4,paste("Control KL replay failure",branch))
  cat("FROZEN_CONTROL_REPLAY_PASS",branch,"Z_MAXABS",err,"MAX_TVD",tvd,"\n")
}
manifest <- fread(file.path(SRC,"461_amended_technical_input_manifest_LOCAL_ONLY.tsv"),sep="\t",quote="")
counts <- c(SRP157974=334L,SRP057793=90L,SRP002628=14L,SRP094795=23L)
assert(nrow(manifest)==461L&&uniqueN(manifest$library_id_LOCAL_ONLY)==461L&&
       all(manifest$input_status=="INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED")&&
       all(manifest$final_clinical_approval=="NO_UNSIGNED")&&
       all(manifest$patient_independence_certified=="NO")&&
       sum(as.integer(manifest$stage2_expected_runs))==483L,"Manifest contract")
summ <- function(x){
  x <- as.numeric(x);assert(all(is.finite(x)),"Nonfinite summary")
  c(median=median(x),p95=unname(quantile(x,.95)),max=max(x))
}
write <- function(x,name)fwrite(x,file.path(FINAL,name),sep="\t")
adequacy_agg <- list();sens_agg <- list();component_agg <- list();sentinels <- list()
for(study in names(counts)){
  m <- manifest[manifest$study==study];N_PAIRS <- as.integer(counts[[study]])
  ids <- sort(as.character(m$library_id_LOCAL_ONLY))
  assert(length(ids)==N_PAIRS,"Study dimensions")
  q <- fread(file.path(SRC,study,"library_expression_QC_LOCAL_ONLY.tsv"),sep="\t",quote="")
  q <- q[match(ids,q$library_id_LOCAL_ONLY)]
  assert(identical(as.character(q$library_id_LOCAL_ONLY),ids)&&all(toupper(as.character(q$UQ_evaluable))=="TRUE"),"QC axis or UQ mismatch")
  store <- list()
  for(branch in BRANCHES){
    X <- read_499(file.path(SRC,study,paste0("499_recount3_",branch,"_LOCAL_ONLY.tsv.gz")),ids)
    Z <- prepare_z(X,branch)
    bundle <- full500_bundle(Z,W)
    H <- fixed500_project_raw_h(Z,W,paste(study,branch,sep="_"))
    f <- score_fractions(H)
    assert(max(abs(rowSums(f)-1))<=SUM_TO_ONE_TOLERANCE,"Fraction closure")
    Vhat <- W%*%t(H);dimnames(Vhat) <- dimnames(bundle$V)
    a <- calculate_reconstruction_adequacy(bundle$V,Vhat,f)
    assert(identical(as.character(a$sample_id),ids),"Adequacy axis")
    metrics <- c("normalized_KL_per_input_mass","reconstruction_RMSE","reconstruction_R2","reconstruction_cosine","normalized_score_entropy")
    assert(all(is.finite(as.matrix(a[,..metrics]))),"Invalid adequacy")
    a[,study:=study];a[,branch:=branch]
    for(cmp in COMPONENTS){
      a[[paste0("fraction_",cmp)]] <- f[,cmp]
      a[[paste0("raw_H_",cmp)]] <- H[,cmp]
    }
    write(a,paste0(study,"_",branch,"_projection_QC_LOCAL_ONLY.tsv"))
    for(metric in metrics){
      ss <- summ(a[[metric]])
      adequacy_agg[[length(adequacy_agg)+1L]] <- data.table(study=study,branch=branch,n_libraries=N_PAIRS,metric=metric,median=ss[1],p95=ss[2],maximum=ss[3],minimum=min(a[[metric]]),constant_W_rows=length(bundle$zero_variance_W_features))
    }
    for(cmp in COMPONENTS)component_agg[[length(component_agg)+1L]] <- data.table(study=study,branch=branch,component=cmp,n_libraries=N_PAIRS,median_fraction=median(f[,cmp]),dominant_count=sum(dominant_component(f)==cmp))
    store[[branch]] <- list(Z=Z,H=H,f=f,adequacy=a,constant_W_rows=bundle$zero_variance_W_features)
  }
  for(alt in c("TMM_CPM","UQ_CPM")){
    tag <- paste0("CPM_vs_",alt)
    s <- normalization_sensitivity(store$CPM$f,store[[alt]]$f,tag)
    write(s$sample,paste0(study,"_",tag,"_LOCAL_ONLY.tsv"))
    v <- s$sample$total_variation_distance;ss <- summ(v)
    sens_agg[[length(sens_agg)+1L]] <- data.table(study=study,comparison=tag,n_libraries=N_PAIRS,median_TVD=ss[1],p95_TVD=ss[2],max_TVD=ss[3],dominant_agreement=mean(s$sample$dominant_component_concordant))
  }
  # Deterministic first 12 IDs; compare to original full-study result as well as repeated reversed fit.
  z0 <- store$CPM$Z[,ids[seq_len(12L)],drop=FALSE]
  f0 <- store$CPM$f[colnames(z0),,drop=FALSE]
  f1 <- score_fractions(fixed500_project_raw_h(z0,W,paste0(study,"_REPEAT")))
  f2 <- score_fractions(fixed500_project_raw_h(z0[,rev(colnames(z0)),drop=FALSE],W,paste0(study,"_REVERSE")))
  f2 <- f2[rownames(f1),,drop=FALSE]
  d <- rowSums(abs(f1-f2))/2
  delta <- max(abs(f0-f1))
  assert(delta<=1e-12&&median(d)<=1e-5&&max(d)<=1e-3,"Repeat/order sentinel")
  sentinels[[study]] <- data.table(study=study,n=12L,repeat_max_abs_fraction=delta,reverse_median_TVD=median(d),reverse_max_TVD=max(d))
  saveRDS(store,file.path(FINAL,paste0(study,"_fixed500_checkpoint_LOCAL_ONLY.rds")),compress="gzip")
  cat("STUDY_FIXED500_COMPLETE",study,N_PAIRS,"\n")
}
write(rbindlist(adequacy_agg),"aggregate_reconstruction_adequacy.tsv")
write(rbindlist(sens_agg),"aggregate_normalization_sensitivity.tsv")
write(rbindlist(component_agg),"aggregate_component_distribution_DESCRIPTIVE.tsv")
write(rbindlist(sentinels),"aggregate_determinism.tsv")
# Descriptive original ARCHS4 primary adequacy reference by actual frozen cancer class.
ref <- fread(find_frozen(S10,"19_primary_CPM_fixed_W_model_adequacy.tsv.gz"),sep="\t",quote="")
canon <- fread(find_frozen(S9,"01_canonical_library_manifest.tsv.gz"),sep="\t",quote="")
assert(nrow(ref)==9465L&&nrow(canon)==9465L&&uniqueN(ref$sample_id)==9465L&&uniqueN(canon$sample_id)==9465L,"Original atlas axes")
ii <- match(ref$sample_id,canon$sample_id);assert(!anyNA(ii),"Original atlas identity")
ref[,cancer:=as.character(canon$analysis_cancer_class[ii])]
metrics <- c("normalized_KL_per_input_mass","reconstruction_RMSE","reconstruction_R2","reconstruction_cosine")
rr <- rbindlist(lapply(metrics,function(metric){
  assert(metric %in% names(ref),"Reference metric missing")
  ref[, {v<-get(metric);ok<-is.finite(v);list(metric_name=metric,n_total=.N,n_finite=sum(ok),median=if(any(ok))median(v[ok]) else NA_real_,p95=if(any(ok))unname(quantile(v[ok],.95)) else NA_real_)},by=cancer]
}))
write(rr,"aggregate_ARCHS4_original_CPM_reference_DESCRIPTIVE.tsv")
cat("V8_JOINT461_FIXED500_COMPATIBILITY_COMPLETE_PENDING_REVIEW\n")
cat("UNMATCHED_STUDIES_NOT_ASSAY_EQUIVALENCE_NO_CLINICAL_ADMISSIONS\n")
