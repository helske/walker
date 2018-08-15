.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  if(length(stanmodels) > 0) ## for testing purposes
   for (m in modules) {
     loadModule(m, what = TRUE)
   }
}
