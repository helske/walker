.onAttach <- function(libname, pkgname) {
  note <- "Note: Since walker version 1.0.0, the prior distribution for the standard deviation parameters is Gamma(shape, rate)."
  packageStartupMessage(paste(strwrap(note), collapse = "\n"))
}