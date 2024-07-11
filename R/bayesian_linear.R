#' Bayesian Linear Regression
#'
#' @param X A matrix of predictors.
#' @param y A response variable.
#' @param priors A list containing priors for the coefficients and variance (mean and precision matrix).
#' @param burn_in The number of iterations to discard at the beginning of the MCMC chain.
#' @param n_samples The number of samples to draw from the posterior distribution.
#' @param verbose Print iteration details if TRUE.
#' @return A list containing the estimated coefficients, credible intervals, and other statistics.
#' @export
bayesian_linear_regression <- function(X, y, priors, burn_in = 1000, n_samples = 5000, verbose = FALSE) {
  # Add intercept term
  X <- cbind(1, X)
  
  # Number of observations and predictors
  N <- nrow(X)
  p <- ncol(X)
  
  # Prior parameters
  beta_prior_mean <- priors$beta_prior_mean
  beta_prior_precision <- priors$beta_prior_precision
  sigma2_prior_shape <- priors$sigma2_prior_shape
  sigma2_prior_rate <- priors$sigma2_prior_rate
  
  # Initialize storage for MCMC samples
  beta_samples <- matrix(0, n_samples, p)
  sigma2_samples <- numeric(n_samples)
  
  # Initial values
  beta <- rep(0, p)
  sigma2 <- 1
  
  # Gibbs sampling
  for (iter in 1:(burn_in + n_samples)) {
    # Sample from the conditional posterior of beta
    beta_cov <- solve(t(X) %*% X / sigma2 + beta_prior_precision)
    beta_mean <- beta_cov %*% (t(X) %*% y / sigma2 + beta_prior_precision %*% beta_prior_mean)
    beta <- MASS::mvrnorm(1, beta_mean, beta_cov)
    
    # Sample from the conditional posterior of sigma2
    shape <- sigma2_prior_shape + N / 2
    rate <- sigma2_prior_rate + sum((y - X %*% beta)^2) / 2
    sigma2 <- 1 / rgamma(1, shape, rate)
    
    # Store samples after burn-in
    if (iter > burn_in) {
      beta_samples[iter - burn_in, ] <- beta
      sigma2_samples[iter - burn_in] <- sigma2
    }
    
    if (verbose && iter %% 1000 == 0) {
      cat("Iteration:", iter, "\n")
    }
  }
  
  # Calculate credible intervals for the coefficients
  beta_mean <- colMeans(beta_samples)
  beta_sd <- apply(beta_samples, 2, sd)
  credible_intervals <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))
  
  # Create a table for coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(X)[-1]),
    Estimate = beta_mean,
    `Std. Error` = beta_sd,
    `2.5%` = credible_intervals[1, ],
    `97.5%` = credible_intervals[2, ]
  )
  rownames(coefficients_table) <- NULL
  
  # Create a table for sigma2
  sigma2_mean <- mean(sigma2_samples)
  sigma2_sd <- sd(sigma2_samples)
  sigma2_ci <- quantile(sigma2_samples, probs = c(0.025, 0.975))
  
  sigma2_table <- data.frame(
    Parameter = "sigma2",
    Estimate = sigma2_mean,
    `Std. Error` = sigma2_sd,
    `2.5%` = sigma2_ci[1],
    `97.5%` = sigma2_ci[2]
  )
  
  # Return the results as a list
  results <- list(
    coefficients = coefficients_table,
    sigma2 = sigma2_table
  )
  
  return(results)
}
