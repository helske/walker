#' Prediction intervals for walker object
#' 
#' Plots sample quantiles from pre
#' See \code{\link[bayesplot]{ppc_ribbon}} for details.
#' 
#' @importFrom ggplot2 ggplot facet_wrap geom_ribbon geom_line 
#' @importFrom bayesplot color_scheme_get theme_default
#' @param object An output from \code{\link{predict.walker_fit}}.
#' @param level Level for intervals. Default is 0.05, leading to 90\% intervals.
#' @param alpha Transparency level for \code{geom_ribbon}.
#' @export
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
      aes_(~x, ~y, alpha = 0.5), inherit.aes = FALSE) 
  
}