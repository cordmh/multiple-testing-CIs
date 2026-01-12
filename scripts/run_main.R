# === run_main.R ===

# 1. Setup ----
required_packages <- c("ggplot2", "sandwich")
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(5868)
dir.create("outputs", showWarnings = FALSE)

# 2. Source functions ----
source("R/helpers.R")
source("R/clean_data.R")
source("R/algorithm_a.R")
source("R/algorithm_b.R")
source("R/gen_max_ir_dist.R")

# 3. Load input data ----
if (!file.exists("data/funds_data.csv")) {
  stop(paste(
    "Funds dataset not found. Place 'funds_data.csv' in the data/ folder (note", 
    "this is a private data-set, excluded from the public repository).",
    sep = " "
  ))
}
if (!file.exists("data/factors_data.csv")) {
  stop("Factors dataset not found. Place 'factors_data.csv' in the data/ folder.")
}
funds_data <- read.csv("data/funds_data.csv", na.strings = "")
factors_data <- read.csv("data/factors_data.csv")

# 4. Data cleaning ----
cleaned_data <- clean_data(funds_data, factors_data)

# 5. Run main algorithms ----
null_alpha_dist <- run_algorithm_a(cleaned_data)
confidence_intervals <- run_algorithm_b(cleaned_data)
max_ir_distribution <- gen_max_ir_dist(cleaned_data)

# 6. Done ----
cat("Pipeline completed successfully.\n")

# === End of run_main.R ===