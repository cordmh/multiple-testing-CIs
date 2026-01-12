#' Apply Algorithm B from Cord (2025)
#'
#' This function applies Algorithm B (Cord 2025) given a matrix of fund returns and a matrix of
#' factor returns. Or more generally, given a matrix of response variables and a matrix of
#' explanatory variables where each column of the former is regressed on the latter.
#'
#' @data A list of two elements containing a matrix of fund returns and a matrix of factor
#'   returns respectively.
#' @param perf_measure A string specifying the performance measure used to evaluate the funds,
#'   either `alpha` or `ir`. Default is `alpha`.
#' @param sample_rank An integer specifying the rank of the fund to evaluate, or a string equal
#'   to `all` to return results for all ranks. Default is `all`.
#' @param B An integer specifying the number of outer bootstrap simulations for each fund.
#'   Default is 5000.
#' @param S An integer specifying the number of inner bootstrap simulations for each fund.
#'   Default is 500.
#' @param conf_level A numeric value specifying the desired confidence level. Default is 0.9.
#' @param plot A logical value specifying whether to plot the simulated t-distribution
#' @return ci A data frame containing the confidence interval/s and its components.
#' @examples
#' run_algorithm_b(data)
#'
run_algorithm_b <- function(data, perf_measure = "alpha", sample_rank = "all", B = 2500, S = 250, conf_level = 0.9, plot = FALSE) {
  # Initialise data structures:
  start_time <- Sys.time()
  funds <- data[[1]]
  factors <- data[[2]]
  fundind <- nonmissing_indices(funds) # as in algorithm A
  N <- ncol(funds)
  bs_ind_array <- vector("list", length = N) # matrix of indices
  alpha <- vector(length = N)
  alpha_star <- matrix(nrow = N, ncol = B)
  alpha_star_star <- array(0, dim = c(N,B,S))
  se_alpha <- vector(length = N)
  ci <- data.frame(matrix(nrow = N, ncol = 7))
  colnames(ci) <- c("sample_rank", "estimate", "se", "t_lower", "t_upper", "ci_lower", "ci_upper")

  # Perform the iterated bootstrap procedure:
  for (i in 1:N) {
    # step 1:
    funds_i <- funds[fundind[1,i]:fundind[2,i],][,i]
    factors_i <- factors[fundind[1,i]:fundind[2,i],]
    Ti <- length(funds_i) # no. returns for ith fund
    reg <- lm(funds_i ~ factors_i)
    alpha[i] <- reg$coef[[1]]
    se_alpha[i] <- as.numeric(sqrt(diag(kernHAC(reg))[1]))
    # ir[i] = alpha[i]/summary(reg)$sigma
    # step 2:
    bs_ind_array[[i]] <- replicate(B, sample(1:Ti, replace = TRUE), simplify = "matrix")
    for (b in 1:B) {
      # step 3:
      ind_star <- bs_ind_array[[i]][,b]
      factors_star <- factors_i[ind_star,]
      return_star <- funds_i[ind_star]
      # step 4:
      alpha_star[i,b] <- lm(return_star ~ factors_star)$coef[[1]]
      for (s in 1:S) {
        # step 5:
        ind_star_star <- sample(ind_star, replace = TRUE)
        # step 6:
        return_star_star <- funds_i[ind_star_star]
        factors_star_star <- factors_i[ind_star_star,]
        # step 7:
        alpha_star_star[i,b,s] <- lm(return_star_star ~ factors_star_star)$coef[[1]]
      }
    }
    cat(paste0("  Fund ", i, " done\n"))
  }

  # Process bootstrap results:
  if (is.numeric(sample_rank)) {
    grid <- as.integer(sample_rank)
  } else {
    grid <- 1:N
  }
  for (k in grid) {
    # step 8:
    k_alpha <- kth_largest(alpha, k)
    std_err_k <- se_alpha[which_kth_largest(alpha, k)]
    # step 9:
    k_alpha_star <- apply(alpha_star, 2, kth_largest, k)
    # step 10:
    k_alpha_bc <- mean(k_alpha_star)
    # step 11:
    k_alpha_star_star <- apply(alpha_star_star, c(2,3), kth_largest, k)
    # step 12:
    se_k_alpha_star <- apply(k_alpha_star_star, 1, sd)
    # step 13:
    t <- (k_alpha_star - k_alpha) / se_k_alpha_star
    # step 14:
    lower_ptile <- (1 - conf_level) / 2
    upper_ptile <- (conf_level + 1) / 2
    t_lower <- as.numeric(quantile(t, lower_ptile))
    t_upper <- as.numeric(quantile(t, upper_ptile))
    # step 15:
    ci[k, 1] <- k
    ci[k, 2] <- k_alpha_bc
    ci[k, 3] <- std_err_k
    ci[k, 4] <- t_lower
    ci[k, 5] <- t_upper
    ci[k, 6] <- k_alpha_bc - t_upper * std_err_k
    ci[k, 7] <- k_alpha_bc - t_lower * std_err_k
    # Plot bootstrap-t distribution/s:
    if (plot) {
      segments_df <- data.frame(
        x = c(x = t_lower, x = t_upper),
        xend = c(t_lower, t_upper),
        y = c(0, 0),
        yend = c(Inf, Inf)
      )
      temp <- ifelse(perf_measure == "alpha", "Alpha ", "Information ratio ")
      plot_title <- paste0(temp, "bootstrap-t distribution (from iterated bootstrap procedure)")
      ggplot(as.data.frame(t), aes(x = t)) +
        geom_histogram(aes(y = after_stat(density))) +
        labs(x = "Bootstrap t-statistic", y = "Density", title = plot_title) +
        geom_segment(data = segments_df, aes(x = x, xend = xend, y = y, yend = yend))
      ggsave(paste0("outputs/", perf_measure, "_", k, "_t_dist.png"), width = 20, height = 10, units = "cm", dpi = 300)
    }
  }

  # Print confidence interval:
  if (is.numeric(sample_rank)) {
    temp1 <- as.character(round(k_alpha_bc, 2))
    temp2 <- as.character(round(t_lower, 2))
    temp3 <- as.character(round(std_err_k, 2))
    temp4 <- as.character(round(t_upper, 2))
    temp5 <- paste0("\n  The ", conf_level * 100, "% CI is: [", temp1, "-", temp4, "*", temp3, ",", temp1, "-", temp2, "*", temp3, "]")
    temp6 <- gsub("--", "+", temp5)
    temp7 <- paste0(" = [", round(ci[k,6], 3), ",", round(ci[k,7], 3), "].\n\n")
    cat(paste0(temp6, temp7)) # print confidence interval and its components
  }
  mins_elapsed <- round((as.numeric(difftime(Sys.time(), start_time, units = "mins"))), 0)
  temp <- paste0("\n  It took ", mins_elapsed, " minutes to run the algorithm.\n\n")
  cat(sub("1 minutes", "1 minute", temp))
  cat(paste0("  Enjoy your ", get_time_period(), "!\n\n"))


  # Return confidence interval/s:
  if (is.numeric(sample_rank)) {
    out <- ci[sample_rank, ]
  } else {
    out <- ci
  }
  return(out)
}
