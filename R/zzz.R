## Unlooad compiled code when unloading FRK
.onUnload <- function (libpath) {
  library.dynam.unload("FRK", libpath)
}
