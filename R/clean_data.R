# clean_data.R
# Function to load, clean, and preprocess raw input data.
clean_data <- function(funds_data, factors_data) {

  ## Read in and clean factor returns data:
  temp <- tempfile()
  download.file("http://homepage.sns.it/marmi/files/Australia.zip", temp)
  factors_aus_4 <- read.table(unz(temp, "Australia_Factors.txt"), fill = TRUE)
  unlink(temp)
  factors_aus_2 <- factors_data[12:279, ][, c(2:6)]
  rf <- as.numeric(as.matrix(factors_aus_4[36:303, ][, 3]))
  factors_aus_1 <- as.numeric(as.matrix(factors_aus_4[36:303, ][, c(7)]))
  factors_aus_5 <- as.numeric(as.matrix(factors_aus_2)) * 100
  factors <- c(factors_aus_5, factors_aus_1)
  dim(factors) <- c(268, 6)
  colnames(factors) <- c("RMRF", "SMB", "HML", "RMW", "CMA", "MOM")
  factors <- factors[, !colnames(factors) %in% c("RMW", "CMA")] # use Fama-French (1993) and Carhart (1997) factors only.

  ## Read in and clean fund returns data:
  funds <- funds_data
  # Remove irrelevant rows:
  time <- funds[1]
  funds <- data.frame(funds)
  funds <- funds[suppressWarnings(grepl("^a$", unlist(funds[1, ])))]
  funds <- cbind(time, funds)
  funds <- funds[47:nrow(funds), ][, -c(1)] # Drop irrelevant rows and time periods that extend past factor data
  rownames(funds) <- c()
  funds = as.data.frame(lapply(funds, as.numeric))
  funds[,34] = NULL # noncontinuous returns (possible endogeneity) so manually removed this fund
  # Remove funds with <36 months of returns:
  funds <- funds[, colSums(is.na(funds)) < nrow(funds) - 36]
  # Convert to matrix:
  funds2 <- matrix(as.numeric(unlist(funds)), nrow = nrow(funds))
  funds <- data.frame(funds[27:294, ]) # deleting some returns to ensure months consistent with factors data.
  # Subtract risk-free rate from all fund returns:
  for (i in 1:ncol(funds)) {
    for (m in 1:nrow(funds)) {
      if (!is.na(funds[m, i])) {
        funds[m, i] <- funds[m, i] - rf[m]
      }
    }
  }
  funds <- as.matrix(funds)
  row.names(funds) <- NULL
  # Remove funds that are effectively passive:
  cor_mat <- matrix(nrow = ncol(funds), ncol = ncol(factors))
  for (i in 1:ncol(funds)) {
    ind <- which(!is.na(funds[, i]))
    fund_i <- funds[ind, i]
    for (j in 1:ncol(factors)) {
      cor_mat[i, j] <- cor(funds[ind, i], factors[ind, j])
    }
  }
  funds_to_delete <- which(apply(cor_mat, 1, value_above_threshold, 0.9))
  funds <- funds[, -funds_to_delete]
  # Remove funds that are highly correlated with each-other (e.g. retail/ws, herding):
  nonmiss <- colSums(!is.na(funds))
  cor_mat <- cor(funds, use = "pairwise.complete.obs")
  high_cor_pairs <- which(abs(cor_mat) > 0.95 & upper.tri(cor_mat), arr.ind = TRUE)
  drop_cols <- integer(0)
  for (k in seq_len(nrow(high_cor_pairs))) {
    i <- high_cor_pairs[k, 1]
    j <- high_cor_pairs[k, 2]
    if (i %in% drop_cols || j %in% drop_cols) next
    # remove fund with lower number of observations:
    if (nonmiss[i] < nonmiss[j]) {
      drop_cols <- c(drop_cols, i)
    } else {
      drop_cols <- c(drop_cols, j)
    }
  }
  funds <- funds[, -drop_cols, drop = FALSE]
  cat(paste0("\nData cleaning completed!\n  No. funds: ", ncol(funds), ",\n  No. months: ", nrow(funds), ", \n  No. fund return-month observations: ", sum(!is.na(funds)), ".\n\n"))

  ## Return cleaned data:
  return(list(funds, factors))
}
