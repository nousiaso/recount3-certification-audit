#!/usr/bin/env Rscript
# V7L: read-only direct implementation inspection of the lower-level NMF solver.
# Does not run DeconvoluteSE or .nmf.predict; no study expression is read.
suppressPackageStartupMessages(library(SpatialEcoTyper))
stopifnot(identical(as.character(packageVersion("SpatialEcoTyper")), "1.0.3"))
ns <- asNamespace("SpatialEcoTyper")
cat("V7L_INSTALLED_NMF_PREDICTOR_IMPLEMENTATION_V1\n")
cat("VERSION SpatialEcoTyper ",as.character(packageVersion("SpatialEcoTyper")),"\n",sep="")
cat("R_VERSION ",as.character(getRversion()),"\n",sep="")
cat("NAMESPACE ",environmentName(ns),"\n",sep="")
cat("DIRECT_NMF_PREDICT_EXISTS ",exists(".nmf.predict",envir=ns,inherits=FALSE),"\n",sep="")
if (!exists(".nmf.predict",envir=ns,inherits=FALSE)) {
  cat("LOW_LEVEL_SOLVER_NOT_EXPORTED_OR_NOT_NAMESPACE_LOCAL\n")
  quit(save="no",status=2L)
}
f <- get(".nmf.predict",envir=ns,inherits=FALSE)
stopifnot(is.function(f))
cat("SOLVER_ENVIRONMENT ",environmentName(environment(f)),"\n",sep="")
cat("SOLVER_FORMALS\n")
print(formals(f))
emit <- function(fun,name,limit=480L) {
  code <- deparse(fun,width.cutoff=200L)
  cat("FUNCTION_BEGIN ",name," DEPARSED_LINES ",length(code),"\n",sep="")
  for (i in seq_len(min(length(code),limit))) {
    cat(sprintf("CODE_%s_%03d %s\n",name,i,code[[i]]))
  }
  cat("FUNCTION_END ",name," TRUNCATED ",length(code)>limit,"\n",sep="")
}
emit(f,"LOW_LEVEL_NMF_PREDICT")
if (requireNamespace("codetools",quietly=TRUE)) {
  globals <- codetools::findGlobals(f,merge=FALSE)
  functions <- sort(unique(as.character(globals$functions)))
  variables <- sort(unique(as.character(globals$variables)))
  cat("LOW_LEVEL_REFERENCED_FUNCTION_NAMES ",paste(functions,collapse=";"),"\n",sep="")
  cat("LOW_LEVEL_REFERENCED_VARIABLE_NAMES ",paste(variables,collapse=";"),"\n",sep="")
  local <- functions[vapply(functions,function(nm)
    exists(nm,envir=ns,inherits=FALSE) &&
    is.function(get(nm,envir=ns,inherits=FALSE)),logical(1))]
  cat("LOW_LEVEL_NAMESPACE_FUNCTIONS ",paste(local,collapse=";"),"\n",sep="")
  for (nm in local) {
    if (nm %in% c(".nmf.predict","DeconvoluteSE","NMFpredict")) next
    helper <- get(nm,envir=ns,inherits=FALSE)
    cat("HELPER_FORMALS ",nm," ",paste(names(formals(helper)),collapse=";"),"\n",sep="")
    emit(helper,paste0("HELPER_",gsub("[^[:alnum:]]","_",nm)),limit=90L)
  }
}
cat("PROJECTION_WRAPPER_BRANCHES_FOR_143_AND_9465\n")
for (n in c(143L,9465L)) {
  cat("N_SAMPLES ",n," NSAMPLE_PER_RUN 10000"," SPLIT_BRANCH ",n>10000L,"\n",sep="")
}
cat("NO_EXPRESSION_VALUES_LOADED\nNO_MODEL_PROJECTION\nNO_FROZEN_SOURCE_MODIFIED\nCLINICAL_ADMISSIONS 0\n")
cat("V7L_INSTALLED_NMF_PREDICTOR_IMPLEMENTATION_V1_COMPLETE\n")
