#' Generate the max IR distribution
#'
#' This function approximates the maximum information ratio distribution based on theoretical
#' results derived in Appendix A of Cord (2025).
#'
#' @data A list of two elements containing a matrix of fund returns and a matrix of factor
#'   returns respectively.
#' @param n_iter An integer specifying the number of Monte Carlo simulations for the inverse
#'   CDF method. Default is 50000.
#' @param plot A logical value specifying whether to plot the maximum information ratio
#'   distribution.
#' @return Numeric vector of features of the maximum information ratio distribution.
#'
gen_max_ir_dist <- function(data, n_iter = 50000, plot = TRUE) {
  # initialise data structures:
  funds <- data[[1]]
  factors <- data[[2]]
  fundind <- nonmissing_indices(funds)
  N <- ncol(funds)
  alpha <- vector(length = N) # vector of estimated alphas
  pval <- vector(length = N) # vector of alpha p-values
  ir <- vector(length = N) # vector of information ratios
  sigma <- vector(length = N) # vector of residual standard errors
  Ti_vec <- vector(length = N) # vector of Ti's
  var_alpha <- vector(length = N) # vector of Var(alpha^)s

  # compute relevant parameters for each fund
  for (i in 1:N) {
    funds_i <- funds[fundind[1, i]:fundind[2, i], ][, i]
    factors_i <- factors[fundind[1, i]:fundind[2, i], ]
    Ti <- fundind[2, i] - fundind[1, i] + 1 # no. returns for ith fund
    reg <- lm(funds_i ~ factors_i)
    alpha[i] <- reg$coefficients[["(Intercept)"]]
    ir[i] <- reg$coefficients[["(Intercept)"]] / summary(reg)$sigma
    sigma[i] <- summary(reg)$sigma
    Ti_vec[i] <- length(funds_i)
    temp1 <- cbind(1, reg$model[, c(2)])
    temp2 <- colMeans(temp1)
    temp3 <- solve(t(temp1) %*% temp1)
    var_alpha[i] <- sigma[i] * (1 / Ti + t(temp2) %*% temp3 %*% temp2)
  }

  # use above parameters to form maximum ir distribution:
  ir_expectation <- ir
  ir_variance <- (alpha^2 + var_alpha) / (2 * sigma^2 * (Ti_vec - 5)) + var_alpha / sigma^2
  cdf <- function(x) {
    prod(pnorm((x - ir_expectation) / sqrt(ir_variance)))
  }
  inverse_cdf <- function(p) {
    lower <- min(ir_expectation) - 10 * max(sqrt(ir_variance))
    upper <- max(ir_expectation) + 10 * max(sqrt(ir_variance))
    bounds <- c(lower, upper)
    uniroot(function(x) cdf(x) - p, interval = bounds)$root
  }
  max_ir_dist <- sapply(runif(50000), inverse_cdf)

  if (plot) {
    ggplot(as.data.frame(max_ir_dist), aes(x = max_ir_dist)) +
      geom_histogram(aes(y = after_stat(density))) +
      geom_density() +
      labs(x = "Maximum information ratio", y = "Density", title = "Histogram of maximum information ratio (asymptotic approximation)") +
      xlim(0.6, 1.5)
    ggsave("outputs/max_ir_plot.png", width = 20, height = 10, units = "cm", dpi = 300)
  }

  # quantities of theoretical max ir distribution:
  ptile_2.5th <- unname(quantile(max_ir_dist, 0.025))
  ptile_97.5th <- unname(quantile(max_ir_dist, 0.975))
  expectation <- mean(max_ir_dist)
  variance <- var(max_ir_dist)
  sd <- sqrt(variance)
  bias <- sqrt(sum(ir_variance)) * qnorm(1 - 1 / N)
  mse <- bias^2 + variance
  to_return <- c("ptile_2.5th", "ptile_97.5th", "expectation", "variance", "sd", "bias", "mse")
  return(unlist(mget(to_return)))
}
