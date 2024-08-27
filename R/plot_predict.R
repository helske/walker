#' Prediction intervals for walker object
#' 
#' Plots sample quantiles and posterior means of the predictions 
#' of the `predict.walker_fit` output.
#' 
#' @importFrom ggplot2 ggplot facet_wrap geom_ribbon geom_line 
#' @importFrom bayesplot color_scheme_get theme_default
#' @param object An output from [predict.walker_fit()].
#' @param draw_obs Either `"response"`, `"mean"`, or `"none"`, 
#' where `"mean"` is response variable divided by number of trials or exposures 
#' in case of binomial/poisson models. 
#' @param level Level for intervals. Default is 0.05, leading to 90% intervals.
#' @param alpha Transparency level for [ggplot2::geom_ribbon()].
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
#' 
#' # note very small number of iterations for the CRAN checks!
#' rw2_fit <- walker(y ~ 1 + 
#'                     rw2(~ -1 + x,
#'                         beta = c(0, 10), 
#'                         nu = c(0, 10)),
#'                   beta = c(0, 10), data = data_old,
#'                   iter = 300, chains = 1, init = 0, refresh = 0)
#' 
#' pred <- predict(rw2_fit, newdata = data.frame(x=x[(n-9):n]))
#' data_new <- data.frame(t = (n-9):n, y = y[(n-9):n])
#' plot_predict(pred) + 
#'   ggplot2::geom_line(data = data_new, ggplot2:: aes(t, y), 
#'   linetype = "dashed", colour = "red", inherit.aes = FALSE)
#'
plot_predict <- function(object, draw_obs = NULL, level = 0.05, alpha = 0.33){
  
  
  if (missing(draw_obs)) {
    if(attr(object, "type") == "link") {
      draw_obs <- "none"
    } else {
      if(attr(object, "type") == "mean") {
        draw_obs <- "mean"
      } else draw_obs <- "response" 
    }
  } else {
    draw_obs <- match.arg(draw_obs, c("mean", "response", "none"))
  }
  pred_data <- as.data.frame(as.table(object$y_new))  
  pred_data$time <- as.numeric(levels(pred_data$time))[pred_data$time]
  names(pred_data)[3] <- "value"
  quantiles <- summarise(group_by(pred_data, time), 
    lwr = quantile(.data$value, prob = level), 
    median = quantile(.data$value, prob = 0.5),
    upr = quantile(.data$value, prob = 1 - level))
  
  if (draw_obs != "none") {
    if(draw_obs == "mean" && !is.null(object$u)) {
      obs <- data.frame(y = object$y / object$u, x = time(object$y))
    } else {
      obs <- data.frame(y = object$y, x = time(object$y))
    }
  }
  p <- ggplot(
    data = quantiles,
    mapping = aes(
      x = .data$time,
      y = .data$median,
      ymin = .data$lwr,
      ymax = .data$upr
    )
  )  +
    geom_ribbon(aes_(color = "y_new", fill = "y_new"),
      alpha = alpha, linetype = 0) +
    geom_line(aes_(color = "y_new")) +
    labs(y = NULL) +  theme_default() + theme(legend.position = "none") + 
    scale_fill_manual(
      name = "",
      values = c(y_new = color_scheme_get()[[1]]))
  if(attr(object, "type") != "link") {
    p <- p + geom_line(data = data.frame(y = object$mean, x = time(object$mean)), 
      aes_(~x, ~y, alpha = 1, color = "mean"), inherit.aes = FALSE) + 
      scale_color_manual(
        name = "",
        values = c(y_new = color_scheme_get()[[2]], mean = color_scheme_get()[[4]]))
  } else {
    p <- p +  scale_color_manual(
      name = "",
      values = c(y_new = color_scheme_get()[[2]]))
  }
  if(draw_obs != "none") {
    p <- p + geom_line(data = obs, 
      aes_(~x, ~y, alpha = 1), inherit.aes = FALSE)
  }
  p
}
