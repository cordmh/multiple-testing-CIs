# helpers.R
# Utility functions for data manipulation and transformation.

# Return index of first and last non-missing value of each matrix or data frame column:
nonmissing_indices <- function(df) {
  out <- matrix(nrow = 2, ncol = ncol(df))
  for (i in 1:ncol(df)) {
    NonNAindices <- which(!is.na(df[, i]))
    out[1, i] <- min(NonNAindices)
    out[2, i] <- max(NonNAindices)
  }
  out
}

# Returns the kth-largest element in a vector.
kth_largest <- function(x, rank_ = k) {
  sort(x, decreasing = TRUE)[rank_]
}

# Returns the index of the kth-largest element in a vector.
which_kth_largest <- function(x, rank_ = k) {
  order(x, decreasing = TRUE)[rank_]
}

# Returns TRUE if vector has at least one value exceeding a threshold.
value_above_threshold <- function(vec, threshold) {
  any(vec > threshold)
}

# Returns the vector of Bayesian posterior intercepts (using the g-prior
# and inverse-gamma prior for the parameter vector and error variance
# respectively); where y is the response variable, X is the design matrix,
# S is n_iter; and g, phi_sq, nu are prior distribution hyperparameters.
bayesian_intercepts <- function(y, X, g, phi_sq, nu, S) {
  Hg <- (g / (g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
  SSRg <- t(y) %*% (diag(1, nrow = nrow(X)) - Hg) %*% y
  sigma_sq <- 1 / rgamma(S, (nu + nrow(X)) / 2, (nu * phi_sq + SSRg) / 2)
  coef_var <- g * solve(t(X) %*% X) / (g + 1)
  coef_expectation <- coef_var %*% t(X) %*% y
  E <- matrix(rnorm(S * ncol(X), 0, sqrt(sigma_sq)), S, ncol(X))
  beta <- t(t(E %*% chol(coef_var)) + c(coef_expectation))
  return(beta[, 1])
}

# Map current time to applicable period.
get_time_period <- function() {
  hour <- as.numeric(format(Sys.time(), "%H"))
  cut(hour,
    breaks = c(0, 4, 12, 17, 21, 24),
    labels = c("night", "morning", "afternoon", "evening", "night"),
    right = FALSE
  )
}
