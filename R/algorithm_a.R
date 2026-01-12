#' Apply Algorithm A from Cord (2025)
#'
#' This function applies Algorithm A (Cord 2025) given a matrix of fund returns and a matrix of
#' factor returns. Or more generally, given a matrix of response variables and a matrix of
#' explanatory variables where each column of the former is regressed on the latter.
#'
#' @data A list of two elements containing a matrix of fund returns and a matrix of factor
#'   returns respectively.
#' @param perf_measure A string specifying the performance measure used to evaluate the funds,
#'   either `alpha` or `ir`. Default is `alpha`.
#' @param sample_rank An integer specifying the rank of the fund to evaluate, or a string equal
#'   to `all` to return results for all ranks. Default is `all`.
#' @param B An integer specifying the number of bootstrap simulations. Default is 5000.
#' @param impose_H0 A logical value specifying whether the joint null should be imposed (True if
#'   a hypothesis test is required; False if a confidence interval is required, in which case
#'   the function acts solely as a helper function for algorithm_B). Default is TRUE.
#' @param bayesian A logical value specifying whether to use Bayesian regression (simulating from
#'   the posterior distribution) rather than bootstrap simulation. Default is FALSE.
#' @param plot A logical value specifying whether the results of the hypothesis test (simulated
#'   values under H0 vs. estimated value) should be plotted. Only valid if `impose_H0` is TRUE
#'   and `sample_rank` is an integer (otherwise, will be ignored).
#' @param plot_ranks A numeric vector specifying the ranks of the funds to plot when `bayesian`
#'   is TRUE and `plot` is TRUE. Default is NULL.
#' @return out A vector of simulated alphas (if sample_rank is an integer) or a matrix of simulated alphas
#'   (if sample_rank is `all`).
#' @examples
#' run_algorithm_a(data)
#'
run_algorithm_a <- function(data, perf_measure = "alpha", sample_rank = "all", B = 5000, impose_H0 = TRUE, bayesian = FALSE, plot = TRUE, plot_ranks = NULL) {
  # initialise data structures:
  start_time <- Sys.time()
  funds <- data[[1]]
  factors <- data[[2]]
  fundind <- nonmissing_indices(funds) # fundind[1,i] and fundind[2,i] is T_{i,1} and T{i,T_{i}} as defined in report.
  N <- ncol(funds)
  alpha <- vector(length = N) # vector of estimated alphas
  pval <- vector(length = N) # vector of alpha p-values
  ir <- vector(length = N) # vector of information ratios
  sigma <- vector(length = N) # vector of residual standard errors
  Ti_vec <- vector(length = N) # vector of Ti's
  var_alpha <- vector(length = N) # vector of Var(alpha^)s
  bs_ind_array <- vector("list", length = N) # matrix of indices
  bs <- matrix(nrow = B, ncol = N) # bootstrapped values of the chosen performance measure

  # steps 1-5:
  for (i in 1:N) {
    # initialise:
    funds_i <- funds[fundind[1, i]:fundind[2, i], ][, i]
    factors_i <- factors[fundind[1, i]:fundind[2, i], ]
    Ti <- length(funds_i) # no. returns for ith fund
    # step 1:
    reg <- lm(funds_i ~ factors_i)
    sigma[i] <- summary(reg)$sigma
    if (bayesian) {
      bs[, i] <- bayesian_intercepts(funds_i, cbind(1, factors_i), Ti, sigma[i]^2, 1, B)
      alpha[i] <- mean(bs[, i]) # posterior mean
      if (impose_H0) bs[, i] <- bs[, i] - mean(bs[, i])
      if (perf_measure != "alpha") {
        stop("For Bayesian regression, `perf_measure` must be `alpha`.")
      }
    } else {
      alpha[i] <- reg$coefficients[["(Intercept)"]]
      ir[i] <- reg$coefficients[["(Intercept)"]] / summary(reg)$sigma
      bs_ind_array[[i]] <- replicate(B, sample(1:Ti, replace = TRUE), simplify = "matrix")
      for (b in 1:B) {
        # step 2:
        ind <- bs_ind_array[[i]][, b]
        bs_factors <- factors_i[ind, ]
        bs_return <- funds_i[ind]
        # step 3:
        if (impose_H0) bs_return <- bs_return - alpha[i]
        bs_reg <- lm(bs_return ~ bs_factors)
        # step 4:
        if (perf_measure == "alpha") {
          bs[b, i] <- bs_reg$coefficients[[1]]
        }
        if (perf_measure == "ir") {
          bs[b, i] <- bs_reg$coefficients[[1]] / summary(bs_reg)$sigma
        }
      }
    }
    cat(paste0("  Fund ", i, " done\n"))
  }
  if (perf_measure == "alpha") {
    perf_measures <- alpha
  }
  if (perf_measure == "ir") {
    perf_measures <- ir
  }

  # step 5:
  corr_matrix <- cor(bs)
  ave_abs_corr <- mean(abs(corr_matrix[lower.tri(corr_matrix)]))
  A_star <- t(apply(bs, 1, sort, decreasing = T))

  # step 6:
  if (sample_rank == "all") {
    out <- A_star
  } else {
    out <- A_star[, sample_rank]
  }
  if (!impose_H0 && !plot && is.numeric(sample_rank)) {
    index <<- order(perf_measures, decreasing = TRUE)[sample_rank]
    return(list(out, bs_ind_array))
  }

  # step 7:
  if (impose_H0 && is.numeric(sample_rank)) {
    test_stat <- sort(perf_measures, decreasing = TRUE)[sample_rank]
    pval <- sum(out > test_stat) / B
    str <- if (pval <= 0.05) "skilled" else "not skilled"
  }

  # plot (optional):
  if (plot && impose_H0 && !bayesian && is.numeric(sample_rank)) {
    if (perf_measure == "alpha") {
      x_label_text <- "Alpha (%pm)"
      temp <- paste0("Alpha null distribution: ", sample_rank, "th-ranked fund")
      title_text <- sub("1th", "Highest", temp)
    }
    if (perf_measure == "ir") {
      x_label_text <- "Information ratio"
      temp <- paste0("Information ratio null distribution: ", sample_rank, "th-ranked fund")
      title_text <- sub("1th", "Highest", temp)
    }
    ggplot(as.data.frame(out), aes(x = out)) +
      geom_histogram(aes(y = after_stat(density))) +
      labs(x = x_label_text, y = "Density", title = title_text) +
      geom_segment(aes(x = test_stat, y = 0, xend = test_stat, yend = Inf), linewidth = 1)
    ggsave(paste0("outputs/", perf_measure, "_", sample_rank, "_hypothesis_test.png"), width = 20, height = 10, units = "cm", dpi = 300)
  }
  if (plot && impose_H0 && !is.numeric(sample_rank) && bayesian) {
    message <- " if `bayesian` and `plot` are both TRUE"
    if (is.null(plot_ranks)) stop(paste0("plot_ranks must be specified", message))
    alphas_sorted <- sort(alpha, decreasing = TRUE)
    for (i in seq_along(plot_ranks)) {
      ggplot(as.data.frame(A_star), mapping = aes_string(paste0(x = "V", as.character(i)))) +
        geom_histogram(aes(y = after_stat(density)), bins = 50) +
        geom_segment(aes(x = alphas_sorted[i], y = 0, xend = alphas_sorted[i], yend = Inf), linewidth = 1.2) +
        labs(x = "Alpha", y = "Density")
      ggsave(paste0("outputs/", "bayes", i, ".png"), plot = last_plot(), width = 12, height = 10, units = "cm", dpi = 300)
    }
  }

  # print results:
  if (impose_H0 && is.numeric(sample_rank)) {
    index <- order(perf_measures, decreasing = TRUE)[sample_rank]
    cat(paste0("\n  — Fund ", index, " generated the ", perf_measure, " with sample rank ", sample_rank, "."))
    cat(paste0("\n  — The p-value is ", pval, ", so this fund is ", str, " with 95% confidence (based on its ", perf_measure, ")."))
  }
  time_elapsed <- round(as.numeric(Sys.time() - start_time, units = "secs"), 0)
  cat(paste0("\n  — The algorithm took ", time_elapsed, " seconds to run.\n\n"))


  # return null distribution/s:
  return(out)
}
