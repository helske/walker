context("Test walker")

test_that("arguments work as intended", {
  expect_error(walker("aa"))
  expect_error(walker(rnorm(2) ~ 1:4))
  expect_error(walker(rnorm(10) ~ 1))
  expect_error(walker(y ~ rw1(~1)))
  expect_error(walker(rnorm(10) ~ 1, beta_prior = 0))
})

test_that("we get proper output") {
  y <- 1:10
  expect_error(walker(y ~ rw(~ 1, beta_prior = c(0, 1), sigma_prior = c(0, 1))), 
    sigma_y_prior = c(0, 1))
}
