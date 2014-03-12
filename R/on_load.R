## preliminary code

print.synlik.version <- function()
{ library(help=synlik)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is \"synlik\" ",version,"\n"),sep="")
}


.onAttach <- function(...) { 
  print.synlik.version()
}

.onUnload <- function(libpath) library.dynam.unload("synlik", libpath)

.onLoad <- function(lib,pkg) {
   library.dynam("synlik", pkg, lib)
}

.First.lib <- function(lib, pkg) {
  library.dynam("synlik", pkg, lib)
  print.synlik.version()
}
