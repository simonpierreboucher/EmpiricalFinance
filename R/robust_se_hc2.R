#' Robust Standard Errors (HC2)
#'
#' @param model_results A list containing the results from the `linear_model` function.
#'
#' @return A data frame containing the coefficients, robust standard errors, and p-values.
#' @export
robust_se_hc2 <- function(model_results) {
  if (!is.list(model_results) || !"coefficients" %in% names(model_results)) {
    stop("model_results must be a list containing the results from the `linear_model` function.")
  }
  
  coefficients <- as.vector(model_results$coefficients$Value)
  X <- model_results$model_matrix
  residuals <- model_results$residuals
  
  n <- nrow(X)
  k <- ncol(X)
  
  # Compute hat matrix
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  
  # Compute robust variance-covariance matrix
  u_hat <- residuals / sqrt(1 - diag(H))
  meat <- t(X) %*% diag(u_hat^2) %*% X
  bread <- solve(t(X) %*% X)
  V_hat_HC2 <- bread %*% meat %*% bread
  
  # Compute robust standard errors
  se_robust <- sqrt(diag(V_hat_HC2))
  
  # Compute robust t-values and p-values
  t_values_robust <- coefficients / se_robust
  p_values_robust <- 2 * pt(-abs(t_values_robust), df = n - k)
  
  # Create a data frame for the results
  results <- data.frame(
    Coefficient = model_results$coefficients$Coefficient,
    Value = coefficients,
    `Robust Standard Error` = se_robust,
    `Robust P-value` = p_values_robust
  )
  
  return(results)
}
