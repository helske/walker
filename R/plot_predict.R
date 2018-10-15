#' Prediction intervals for walker object
#' 
#' Plots sample quantiles from pre
#' See \code{\link{ppc_ribbon}} for details.
#' 
#' @importFrom ggplot2 ggplot facet_wrap geom_ribbon geom_line 
#' @importFrom bayesplot color_scheme_get theme_default
#' @param object An output from \code{\link{predict.walker_fit}}.
#' @param level Level for intervals. Default is 0.05, leading to 90\% intervals.
#' @param alpha Transparency level for \code{\link{geom_ribbon}}.
#' @export
#' @examples 
#' set.seed(1)
#' n <- 60
#' slope <- 0.0001 + cumsum(rnorm(n, 0, sd = 0.01))
#' beta <- numeric(n)
#' beta[1] <- 1
#' for(i in 2:n) beta[i] <- beta[i-1] + slope[i-1]
#' ts.plot(beta)                
#' x <- rnorm(n, 1, 0.5)
#' alpha <- 2
#' ts.plot(beta * x)
#' 
#' signal <- alpha + beta * x
#' y <- rnorm(n, signal, 0.25)
#' ts.plot(cbind(signal, y), col = 1:2)
#' data_old <- data.frame(y = y[1:(n-10)], x = x[1:(n-10)])
#' rw2_fit <- walker(y ~ 1 + 
#'                     rw2(~ -1 + x,
#'                         beta_prior = c(0, 10), 
#'                         sigma_prior = c(0, 10), 
#'                         slope_prior = c(0, 10)), 
#'                   sigma_y_prior = c(0, 10), 
#'                   beta_prior = c(0, 10),
#'                   iter = 400, chains = 1, data = data_old)
#' 
#' pred <- predict(rw2_fit, newdata = data.frame(x=x[(n-9):n]))
#' data_new <- data.frame(t = (n-9):n, y = y[(n-9):n])
#' plot_predict(pred) + 
#'   geom_line(data=data_new, aes(t, y), linetype="dashed", colour = "red", inherit.aes = FALSE)
#'
plot_predict <- function(object, level = 0.05, alpha = 0.33){
  
  pred_data <- as.data.frame(as.table(object$y_new))  
  pred_data$time <- as.numeric(levels(pred_data$time))[pred_data$time]
  names(pred_data)[3] <- "value"
  grouped <- group_by(pred_data, time)
  quantiles <- summarise_(grouped, 
    .dots = list(
      lwr = ~quantile(value, prob = level), 
      median = ~quantile(value, prob = 0.5),
      upr = ~quantile(value, prob = 1 - level)))
  
  ggplot(
    data = quantiles,
    mapping = aes_(
      x = ~ time,
      y = ~ median,
      ymin = ~ lwr,
      ymax = ~ upr
    )
  )  +
    geom_ribbon(aes_(color = "y_new", fill = "y_new"),
      alpha = alpha, linetype = 0) +
   geom_line(aes_(color = "y_new")) +
    labs(y = NULL) +  theme_default() + theme(legend.position = "none") + 
    scale_color_manual(
      name = "",
      values = c(y_new = color_scheme_get()[[2]])
    ) +
    scale_fill_manual(
      name = "",
      values = c(y_new = color_scheme_get()[[1]])
    ) + geom_line(data = data.frame(y = object$y, x = time(object$y)), 
      aes_(~x, ~y, alpha = 1), inherit.aes = FALSE) 
  
}