#!/usr/bin/env Rscript
# V7K: read-only inspection of installed SpatialEcoTyper 1.0.3
# DeconvoluteSE and its local function-call dependencies.
# NO model runs, no expression loading, no patient or run IDs.
suppressPackageStartupMessages({
  library(SpatialEcoTyper)
})
stopifnot(as.character(packageVersion("SpatialEcoTyper")) == "1.0.3")
ns <- asNamespace("SpatialEcoTyper")
f <- get("DeconvoluteSE", envir=ns, inherits=FALSE)
stopifnot(is.function(f))

cat("V7K_DECONVOLUTESE_INSTALLED_SOURCE_INSPECTION_V1\n")
cat("PACKAGE SpatialEcoTyper\n")
cat("VERSION ",as.character(packageVersion("SpatialEcoTyper")),"\n",sep="")
cat("PACKAGE_LIBRARY ",find.package("SpatialEcoTyper"),"\n",sep="")
cat("R_VERSION ",as.character(getRversion()),"\n",sep="")
cat("DECONVOLUTESE_FUNCTION_IS_NAMESPACE_FUNCTION ",
    identical(f,SpatialEcoTyper::DeconvoluteSE),"\n",sep="")

emit <- function(fun,label,max_lines=400L) {
  code <- deparse(fun,width.cutoff=200L)
  cat("FUNCTION_BEGIN ",label," DEPARSED_LINES ",length(code),"\n",sep="")
  n <- min(length(code),max_lines)
  for (i in seq_len(n)) {
    cat(sprintf("CODE_%s_%03d %s\n",label,i,code[[i]]))
  }
  cat("FUNCTION_END ",label," DISPLAYED_LINES ",n,
      " TRUNCATED ",n<length(code),"\n",sep="")
}
emit(f,"DeconvoluteSE",max_lines=400L)
cat("FUNCTION_FORMAL_ARGUMENTS\n")
print(formals(f))

if (requireNamespace("codetools",quietly=TRUE)) {
  callnames <- unique(as.character(codetools::findGlobals(f,merge=FALSE)$functions))
  callnames <- sort(callnames)
  cat("DECONVOLUTESE_REFERENCED_FUNCTION_NAMES ",
      paste(callnames,collapse=";"),"\n",sep="")
  local <- callnames[vapply(callnames,function(nm) {
    exists(nm,envir=ns,inherits=FALSE) &&
      is.function(get(nm,envir=ns,inherits=FALSE))
  },logical(1))]
  cat("NAMESPACE_LOCAL_CALLEES ",paste(local,collapse=";"),"\n",sep="")
  # Only show function signatures for local dependencies; no recursive
  # blind traversal, which could produce massive stdout.
  for (name in local) {
    g <- get(name,envir=ns,inherits=FALSE)
    cat("LOCAL_CALLEE_FORMALS ",name," ",
        paste(names(formals(g)),collapse=";"),"\n",sep="")
    lines <- deparse(g,width.cutoff=200L)
    cat("LOCAL_CALLEE_BODY_LINES ",name," ",length(lines),"\n",sep="")
    # Show the complete short helpers; otherwise print first 120 lines
    # to let the next diagnostic target the precise dependency.
    if (length(lines)<=120L) emit(g,name,max_lines=120L)
  }
}
cat("NO_EXPRESSION_VALUES_LOADED\n")
cat("NO_DECONVOLUTESE_EXECUTION\n")
cat("NO_SOURCE_FILES_MODIFIED\n")
cat("NO_CLINICAL_ADMISSIONS\n")
cat("V7K_DECONVOLUTESE_INSTALLED_SOURCE_INSPECTION_V1_COMPLETE\n")
