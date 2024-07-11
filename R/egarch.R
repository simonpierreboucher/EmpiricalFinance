#' Estimer un modèle EGARCH(p, q) en utilisant la méthode du maximum de vraisemblance
#'
#' @param series Un vecteur numérique représentant la série temporelle.
#' @param p Un entier représentant l'ordre du terme ARCH.
#' @param q Un entier représentant l'ordre du terme GARCH.
#' @param start_params Un vecteur de paramètres initiaux pour l'optimisation.
#'
#' @return Une liste contenant les paramètres estimés, les erreurs standard, les p-values, et les statistiques du modèle.
#' @export
estimate_egarch <- function(series, p, q, start_params) {
  log_likelihood <- function(params, series, p, q) {
    omega <- params[1]
    alpha <- params[2:(p + 1)]
    gamma <- params[(p + 2):(2 * p + 1)]
    beta <- params[(2 * p + 2):(2 * p + q + 1)]
    
    n <- length(series)
    innovations <- series - mean(series)
    log_variances <- rep(log(var(innovations)), n)
    
    for (i in (max(p, q) + 1):n) {
      past_innovations <- innovations[(i - p):(i - 1)]
      past_log_variances <- log_variances[(i - q):(i - 1)]
      
      log_variances[i] <- omega + 
        sum(alpha * (abs(past_innovations) - sqrt(2/pi))) +
        sum(gamma * past_innovations) +
        sum(beta * past_log_variances)
    }
    
    conditional_variances <- exp(log_variances)
    log_lik <- -0.5 * sum(log(conditional_variances) + (innovations^2 / conditional_variances))
    return(-log_lik)
  }
  
  opt_result <- optim(par = start_params,
                      fn = log_likelihood,
                      series = series,
                      p = p,
                      q = q,
                      method = "BFGS",
                      hessian = TRUE)
  
  estimates <- opt_result$par
  hessian <- opt_result$hessian
  std_errors <- sqrt(diag(solve(hessian)))
  t_values <- estimates / std_errors
  p_values <- 2 * pt(-abs(t_values), df = length(series) - length(estimates))
  
  result <- list(
    coefficients = estimates,
    std_errors = std_errors,
    p_values = p_values,
    log_likelihood = -opt_result$value,
    convergence = opt_result$convergence
  )
  
  return(result)
}


