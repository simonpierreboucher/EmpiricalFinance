#' Logistic Regression using Maximum Likelihood Estimation
#'
#' @param X A matrix of predictors.
#' @param y A binary response variable.
#' @param tol Tolerance level for convergence.
#' @param max_iter Maximum number of iterations.
#' @param verbose Print iteration details if TRUE.
#' @return A list containing the estimated coefficients, standard errors, t-values, p-values, and global statistics.
#' @export
logistic_regression <- function(X, y, tol = 1e-6, max_iter = 1000, verbose = FALSE) {
  # Add intercept term
  X <- cbind(1, X)
  
  # Initialize coefficients
  beta <- rep(0, ncol(X))
  N <- nrow(X)
  J <- ncol(X)
  iter <- 0
  converged <- FALSE
  
  # Iteratively reweighted least squares (IRLS) algorithm
  while (iter < max_iter && !converged) {
    iter <- iter + 1
    eta <- X %*% beta
    mu <- plogis(eta)
    W <- diag(as.vector(mu * (1 - mu)), nrow = N)
    z <- eta + solve(W) %*% (y - mu)
    
    beta_new <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    
    if (max(abs(beta_new - beta)) < tol) {
      converged <- TRUE
    }
    
    beta <- beta_new
    
    if (verbose) {
      cat("Iteration:", iter, " Beta:", beta, "\n")
    }
  }
  
  # Calculate standard errors
  eta <- X %*% beta
  mu <- plogis(eta)
  W <- diag(as.vector(mu * (1 - mu)), nrow = N)
  cov_beta <- solve(t(X) %*% W %*% X)
  se_beta <- sqrt(diag(cov_beta))
  
  # Calculate t-values and p-values
  t_values <- beta / se_beta
  p_values <- 2 * pt(-abs(t_values), df = N - J)
  
  # Calculate global statistics
  log_likelihood <- sum(y * log(mu) + (1 - y) * log(1 - mu))
  null_deviance <- -2 * sum(y * log(mean(y)) + (1 - y) * log(1 - mean(y)))
  residual_deviance <- -2 * log_likelihood
  r_squared <- 1 - (residual_deviance / null_deviance)
  adj_r_squared <- 1 - ((1 - r_squared) * ((N - 1) / (N - J - 1)))
  
  # Create a table for coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(X)[-1]),
    Estimate = beta,
    `Std. Error` = se_beta,
    `t value` = t_values,
    `Pr(>|t|)` = p_values
  )
  rownames(coefficients_table) <- NULL
  
  # Create a table for global statistics
  global_stats <- data.frame(
    `Log Likelihood` = log_likelihood,
    `Null Deviance` = null_deviance,
    `Residual Deviance` = residual_deviance,
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared
  )
  
  # Return the results as a list
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )
  
  return(results)
}

