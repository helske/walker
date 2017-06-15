context("Test walker")

test_that("arguments work as intended", {
  library(walker)
  expect_error(walker("aa"))
  expect_error(walker(rnorm(2) ~ 1:4))
  expect_error(walker(rnorm(10) ~ 1))
  expect_error(walker(y ~ 1))
  expect_error(walker(rnorm(10) ~ 1, beta_prior = 0))
  
  x <- 1:3
  expect_identical(c(1,1,1,1:3), c(walker(1:3 ~ x, return_x_reg = TRUE)))
})
test_that("stan side works", {
  y <- x <- 1:3
  set.seed(1)
  expect_warning(fit <- walker(y ~ x, beta_prior = cbind(0, c(2, 2)),
    sigma_prior = cbind(0, c(2,2,2)), iter = 10, chains = 1, refresh = 0),NA)
  expect_equivalent(structure(c(0.575776440370937, 0.608739297869922, 0.600646410430753, 
    2.47394156830475, 0.503598122307422), .Dim = 5L, .Dimnames = structure(list(
      iterations = NULL), .Names = "iterations")), 
    extract(fit, pars = "sigma_y")$sigma_y)
  
  set.seed(1)
  expect_warning(fit <- walker(y ~ x, naive = TRUE, beta_prior = cbind(0, c(2, 2)),
    sigma_prior = cbind(0, c(2,2,2)), iter = 10, chains = 1, refresh = 0),NA)
  expect_equivalent(structure(c(1.27535032198269, 1.27535032198269, 1.27535032198269, 
    1.07187289906723, 1.24949754280559), .Dim = 5L, .Dimnames = structure(list(
      iterations = NULL), .Names = "iterations")), 
    extract(fit, pars = "sigma_y")$sigma_y)
})