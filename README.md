# Confidence intervals under multiple testing

This repository contains the R code used to produce the results and figures in the paper:

**"Studentized bootstrap confidence intervals for regression parameters under multiple testing"**  
[SSRN Link](https://ssrn.com/abstract=5978695) 

An earlier version was originally submitted for the major project for MATH5868 Bootstrap and Other Resampling Methods (UNSW, T3/2024). 

---

## Abstract

> Suppose a large number of regression models, each with a different response variable but the same set of explanatory variables, are estimated simultaneously. This paper develops a new method to simulate confidence intervals for regression parameters that adjust for multiple testing and whose coverage achieves second-order accuracy under a reasonable set of assumptions. The method is demonstrated by evaluating the performance of a selection of actively managed Australian equity mutual funds.


---

## Project Structure

multiple-testing-CIs/<br>
├── R/ # All functions (algorithms, data cleaning, helpers)<br>
├── data/ # Data<br>
├── outputs/ # Output files<br>
├── scripts/ # Main execution script<br>
├── README.md<br>


## Data Availability

This project relies on a private, commercially licensed dataset of mutual fund returns that **is not included** in this repository.

- The private dataset (of fund returns) is required to run **all** main algorithms and reproduce the core results.
- The public dataset (of factor returns) is included for demonstration purposes **only** and cannot be used to fully run the analyses or replicate the paper.

If you have access to the private data, place it at:

data/funds.csv

and ensure the path in the code matches.

## How to Run

Since all key functions depend on the private data, **you must have access to it** before running any scripts.

1. Obtain and place the private dataset in data/ as described above.
2. Open `scripts/run_main.R` in R or RStudio.
3. Run the entire script to reproduce the paper’s results.

Without the private data, the code will not execute successfully.

Alternatively, source the main script directly:
```r
source("scripts/run_main.R")


Dependencies

This project uses R version 4.5.1 and the following CRAN packages:

    ggplot2

    sandwich