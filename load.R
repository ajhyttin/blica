#These are some of the parameters needed for operation
#

#loads the whole directory to R
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# since the functions here are not really changed, loading them already here
loud<-function() {
  require('mvtnorm')
  require('gtools')
  require('lpSolve') #this is for mcs calculation!

  #cat('Loading the R-code...\n')
  sourceDir('./',trace=FALSE)

  #cat('Done.\n')
}


