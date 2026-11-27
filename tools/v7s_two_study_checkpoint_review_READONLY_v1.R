#!/usr/bin/env Rscript
# Read-only aggregate comparison of executed lung V7P and prostate V7R.
# No model fitting, private identifiers, writes, or source changes.
options(width=240)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
repo <- file.path(args[1],"RECOUNT3_CERTIFICATION_GIT")
paths <- c(SRP074349=file.path(repo,"private_data/v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1_0_1/projection_checkpoint_v1_0_1_LOCAL_ONLY.rds"),
           SRP157215=file.path(repo,"private_data/v7r_SRP157215_projection_v1/projection_checkpoint_v7r_v1_LOCAL_ONLY.rds"))
expected <- c(SRP074349=143L,SRP157215=22L)
versions <- c(SRP074349="v7p_v1_0_1",SRP157215="v7r_v1")
components <- c("NonSE",paste0("SE",1:9))
branches <- c("CPM","TMM_CPM","UQ_CPM")
q <- function(x,p) unname(quantile(x,p))
safe_cor <- function(a,b,method) {
  if(length(a)<2L || sd(a)==0 || sd(b)==0) return(NA_real_)
  cor(a,b,method=method)
}
margin <- function(f) apply(f,1,function(z) {z<-sort(z,decreasing=TRUE);z[1]-z[2]})
dom <- function(f) colnames(f)[max.col(f,ties.method="first")]
metrics <- c("normalized_KL_per_input_mass","reconstruction_RMSE","reconstruction_R2","reconstruction_cosine")
pooled <- setNames(vector("list",length(branches)),branches)
bias <- list()
for(study in names(paths)) {
  x <- readRDS(paths[[study]])
  n <- expected[[study]]
  stopifnot(identical(x$version,versions[[study]]),setequal(names(x$per_branch),branches),
            length(x$uq_arch)==n,length(x$uq_rec)==n,
            is.logical(x$uq_arch),is.logical(x$uq_rec),!anyNA(x$uq_arch),!anyNA(x$uq_rec))
  if(study=="SRP157215") stopifnot(identical(x$study,study))
  for(branch in branches) {
    b <- x$per_branch[[branch]]
    a <- b$arch$f
    r <- b$recount3$f
    stopifnot(is.matrix(a),is.matrix(r),identical(dim(a),c(n,10L)),
              identical(dimnames(a),dimnames(r)),identical(colnames(a),components),
              !is.null(rownames(a)),!anyDuplicated(rownames(a)),
              all(is.finite(a)),all(is.finite(r)),all(a>=0),all(r>=0),
              max(abs(rowSums(a)-1))<1e-8,max(abs(rowSums(r)-1))<1e-8)
    mask <- if(branch=="UQ_CPM") x$uq_arch & x$uq_rec else rep(TRUE,n)
    stopifnot(sum(mask)>=2L)
    aa <- b$arch$adequacy
    rr <- b$recount3$adequacy
    stopifnot(nrow(aa)==n,nrow(rr)==n,all(metrics %in% names(aa)),all(metrics %in% names(rr)))
    # Original adequacy helper emits rows in the fraction matrix order; verify IDs if present.
    if("sample_id" %in% names(aa)) stopifnot(identical(as.character(aa$sample_id),rownames(a)))
    if("sample_id" %in% names(rr)) stopifnot(identical(as.character(rr$sample_id),rownames(r)))
    a <- a[mask,,drop=FALSE]; r <- r[mask,,drop=FALSE]
    d <- r-a; tvd <- rowSums(abs(d))/2
    da <- dom(a); dr <- dom(r)
    pooled[[branch]][[study]] <- list(tvd=tvd,agree=da==dr)
    cat("\n=== STUDY",study,"BRANCH",branch,"N",nrow(a),"===\n")
    cat("COMPONENT_AGREEMENT\n")
    tab <- data.frame(component=components,
      median_ARCHS4=apply(a,2,median),median_recount3=apply(r,2,median),
      SD_ARCHS4=apply(a,2,sd),
      Pearson=vapply(seq_len(10L),function(j)safe_cor(a[,j],r[,j],"pearson"),numeric(1)),
      Spearman=vapply(seq_len(10L),function(j)safe_cor(a[,j],r[,j],"spearman"),numeric(1)),
      mean_signed_delta=colMeans(d),median_signed_delta=apply(d,2,median),
      median_absolute_delta=apply(abs(d),2,median),
      P95_absolute_delta=apply(abs(d),2,q,p=.95))
    print(tab,row.names=FALSE,digits=7)
    bias[[paste(study,branch)]] <- setNames(tab$median_signed_delta,components)
    cat("DOMINANT_COMPONENT_CONFUSION\n")
    print(table(ARCHS4=factor(da,levels=components),recount3=factor(dr,levels=components)))
    for(flag in c(FALSE,TRUE)) {
      ix <- which((da!=dr)==flag)
      if(length(ix)) cat("DISCORDANT",flag,"N",length(ix),
        "MEDIAN_ARCHS4_MARGIN",median(margin(a)[ix]),
        "MEDIAN_RECOUNT3_MARGIN",median(margin(r)[ix]),
        "MEDIAN_TVD",median(tvd[ix]),"\n")
    }
    cat("PAIRED_RECONSTRUCTION_ADEQUACY\n")
    for(m in metrics) {
      av<-aa[[m]][mask];rv<-rr[[m]][mask]
      stopifnot(all(is.finite(av)),all(is.finite(rv)))
      cat(m,"ARCHS4_MEDIAN",median(av),"RECOUNT3_MEDIAN",median(rv),
          "MEDIAN_PAIRED_DELTA",median(rv-av),"P95_ABS_DELTA",q(abs(rv-av),.95),
          "SPEARMAN",safe_cor(av,rv,"spearman"),"\n")
    }
  }
}
cat("\n=== CROSS_STUDY_SIGNED_BIAS_COMPARISON ===\n")
for(branch in branches) {
  a<-bias[[paste("SRP074349",branch)]];r<-bias[[paste("SRP157215",branch)]]
  print(data.frame(branch=branch,component=components,lung_median_delta=a,
      prostate_median_delta=r,same_sign=sign(a)==sign(r)),row.names=FALSE,digits=7)
}
cat("\n=== POOLED_DESCRIPTIVE_RESULTS_NOT_STUDY_BALANCED ===\n")
for(branch in branches) {
  tvd<-unlist(lapply(pooled[[branch]],function(z)z$tvd),use.names=FALSE)
  agree<-unlist(lapply(pooled[[branch]],function(z)z$agree),use.names=FALSE)
  cat(branch,"N",length(tvd),"MEDIAN_TVD",median(tvd),"P95_TVD",q(tvd,.95),
      "MAX_TVD",max(tvd),"DOMINANT_AGREE_N",sum(agree),"DOMINANT_AGREEMENT",mean(agree),"\n")
}
cat("V7S_TWO_STUDY_CHECKPOINT_REVIEW_COMPLETE\n")
cat("NO_MODEL_FITTING_NO_FILES_MODIFIED_CLINICAL_ADMISSIONS_0\n")
