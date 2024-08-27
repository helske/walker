check_normal <- function(x, name) {
  if(length(x) != 2) 
    stop(paste0("Argument ", name, " should be a vector of length two, ",
    "defining the mean and standard deviation for the Gaussian prior. "))
  if(!(x[2] > 0)) 
    stop(paste0("Prior standard deviation for ", name, " should be positive. "))
}
check_gamma <- function(x, name) {
  if(length(x) != 2) 
    stop(paste0("Argument ", name, " should be a vector of length two, ",
    "defining the shape and rate for the Gamma prior the parameter ", name, ". "))
  if(!all(x > 0)) 
    stop(paste0("Both parameters of the Gamma prior for the parameter ", name, " should be positive. "))
}
