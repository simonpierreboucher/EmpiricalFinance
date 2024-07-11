#' ARMA/ARIMA Model Estimation using Maximum Likelihood
#'
#' @param data A numeric vector representing the time series data.
#' @param p An integer representing the order of the autoregressive part.
#' @param d An integer representing the order of differencing.
#' @param q An integer representing the order of the moving average part.
#'
#' @return A list containing the estimated coefficients, standard errors, t-values, p-values, and log-likelihood.
#' @export
arma_model <- function(data, p = 0, d = 0, q = 0) {
  # Vérifier si les paramètres sont valides
  if (!is.numeric(data) || length(data) <= max(p, q, d)) {
    stop("Invalid data or ARMA/ARIMA orders.")
  }

  # Différencier les données si nécessaire
  if (d > 0) {
    data <- diff(data, differences = d)
  }

  # Fonction de log-vraisemblance
  log_likelihood <- function(params, data, p, q) {
    n <- length(data)
    ar_params <- params[1:p]
    ma_params <- params[(p + 1):(p + q)]
    sigma2 <- params[(p + q + 1)]

    residuals <- rep(0, n)
    log_likelihood <- 0

    for (t in (max(p, q) + 1):n) {
      ar_term <- sum(ar_params * data[(t - 1):(t - p)])
      ma_term <- sum(ma_params * residuals[(t - 1):(t - q)])
      residuals[t] <- data[t] - ar_term + ma_term
      log_likelihood <- log_likelihood - 0.5 * (log(2 * pi) + log(sigma2) + (residuals[t]^2 / sigma2))
    }

    return(-log_likelihood)
  }

  # Estimation des paramètres par optimisation
  init_params <- c(rep(0, p + q), var(data))
  optim_results <- optim(init_params, log_likelihood, data = data, p = p, q = q, method = "BFGS")

  # Extraire les résultats
  ar_params <- optim_results$par[1:p]
  ma_params <- optim_results$par[(p + 1):(p + q)]
  sigma2 <- optim_results$par[(p + q + 1)]
  log_lik <- -optim_results$value

  # Calculer les erreurs standard
  hessian <- optimHess(optim_results$par, log_likelihood, data = data, p = p, q = q)
  cov_matrix <- solve(hessian)
  se_params <- sqrt(diag(cov_matrix))

  # Calculer les p-values
  t_values <- optim_results$par / se_params
  p_values <- 2 * pt(-abs(t_values), df = length(data) - length(optim_results$par))

  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c(paste0("AR", 1:p), paste0("MA", 1:q), "sigma2"),
    Estimate = optim_results$par,
    `Std. Error` = se_params,
    `t value` = t_values,
    `Pr(>|t|)` = p_values
  )
  rownames(coefficients_table) <- NULL

  # Retourner les résultats
  results <- list(
    coefficients = coefficients_table,
    log_likelihood = log_lik,
    residuals = residuals
  )

  return(results)
}
