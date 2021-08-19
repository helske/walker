context("Test walker")

test_that("arguments work as intended", {
  expect_error(walker("aa"))
  expect_error(walker(rnorm(2) ~ 1:4))
  expect_error(walker(rnorm(10) ~ 1))
  expect_error(walker(y ~ rw1(~1)))
  expect_error(walker(rnorm(10) ~ 1, beta = 0))
})

test_that("we get proper output", {
  y <- 1:10
  expect_error(fit <- walker(y ~ -1 + 
      rw1(~ 1, beta = c(0, 1), sigma = c(2, 1)), 
    sigma_y = c(2, 1), iter = 10, refresh=0), NA)
  expect_s4_class(fit$stanfit, "stanfit")
})

test_that("we get proper output from glm", {
  y <- 1:10
  expect_error(fit <- walker_glm(y ~ -1 + 
      rw1(~ 1, beta = c(0, 1), sigma = c(2, 1)), 
    distribution = "poisson", iter = 10, refresh=0), NA)
  expect_s4_class(fit$stanfit, "stanfit")
})
