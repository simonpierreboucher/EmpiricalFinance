#' Weighted Least Squares (WLS) Estimation
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#' @param weights A vector of weights to be used in the estimation.
#'
#' @return A list containing the estimated coefficients, standard errors, p-values, adjusted R-squared, F-statistic, and p-value of the F-test.
#' @export
weighted_least_squares <- function(data, weights) {
  # Check if the data frame contains at least two columns
  if (ncol(data) < 2) {
    stop("The data frame must contain at least two columns: one dependent variable and at least one independent variable.")
  }

  # Check if the length of weights matches the number of observations
  if (length(weights) != nrow(data)) {
    stop("The length of the weights vector must match the number of observations in the data.")
  }

  # Extract dependent and independent variables
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Add a column of 1s for the intercept

  # Apply weights
  W <- diag(weights)

  # Estimate coefficients using WLS
  beta_wls <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  residuals <- y - X %*% beta_wls

  # Compute standard errors
  Xt_W_X_inv <- solve(t(X) %*% W %*% X)
  se_beta_wls <- sqrt(diag(Xt_W_X_inv))

  # Compute p-values
  t_values <- beta_wls / se_beta_wls
  p_values <- 2 * pt(-abs(t_values), df = nrow(X) - ncol(X))

  # Compute R-squared and adjusted R-squared
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))

  # Compute F-statistic
  MS_reg <- (SS_tot - SS_res) / (ncol(X) - 1)
  MS_res <- SS_res / (nrow(X) - ncol(X))
  f_statistic <- MS_reg / MS_res
  f_pvalue <- pf(f_statistic, ncol(X) - 1, nrow(X) - ncol(X), lower.tail = FALSE)

  # Create a data frame for the coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(data)[-1]),
    Value = beta_wls,
    `Standard Error` = se_beta_wls,
    `P-value` = p_values
  )
  rownames(coefficients_table) <- NULL

  # Create a data frame for the global statistics
  global_stats <- data.frame(
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared,
    `F-statistic` = f_statistic,
    `F-test P-value` = f_pvalue
  )

  # Create a list to store the results
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )

  return(results)
}
