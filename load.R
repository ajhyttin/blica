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
  #require('tmvtnorm')
  require('gtools')
  #cat('Loading packages...\n')
  #require('combinat')
  #require('ggm')
  #require('gtools')
  require('lpSolve') #this is for mcs calculation!

  #This should start a ps-file viewer
  #Replace this if not working
  #if ( .Platform$OS.type == 'unix' ) {
  #  psviewer<<-'gv'
  #} else if ( .Platform$OS.type == 'windows' ) {
  #  psviewer<<-'C:/\"Program Files\"/Ghostgum/gsview/gsview32 -e'
  #} else if ( .Platform$OS.type == 'mac' ) {
  #  psviewer<<-'open -a Preview' 
  #} else {
  #  psviewer<<-NA
  #}
  #cat('PS viewer:',psviewer,'\n')

  #here checking the sachs directory
  #sachs_dir<<- './../simulations/sachsdata/data/'
  #cat('Sachs data dir:',sachs_dir)

  #if ( file.access(paste(sep='',sachs_dir,'/experiment1.csv'),0) == 0 ) {
  #  cat(' is ok.\n')
  #} else {
  #  cat(' is not OK. Sachs data tests will not work.\n')
  #}

  #cat('Loading the R-code...\n')
  sourceDir('./',trace=FALSE)

  #cat('Done.\n')
}


