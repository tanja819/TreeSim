.onLoad <- function(lib, pkg)
{
  library.dynam("TreeSim", pkg, lib)
}

.onUnload <- function(lib)
{
  library.dynam.unload("TreeSim", lib)
}