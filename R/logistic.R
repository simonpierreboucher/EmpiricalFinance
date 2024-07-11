#' Régression logistique par maximum de vraisemblance
#'
#' @param data Un data frame où la première colonne est la variable dépendante (binaire) et les colonnes restantes sont les variables indépendantes.
#' @param start_params Un vecteur de paramètres initiaux pour l'optimisation.
#'
#' @return Une liste contenant les paramètres estimés, les erreurs standard, les p-values, et les statistiques du modèle.
#' @export
logistic_regression <- function(data, start_params = NULL) {
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Fonction de log-vraisemblance
  log_likelihood <- function(params, X, y) {
    beta <- params
    LP <- X %*% beta
    mu <- plogis(LP)
    L <- dbinom(y, size = 1, prob = mu, log = TRUE)
    -sum(L)
  }
  
  # Paramètres initiaux si non fournis
  if (is.null(start_params)) {
    start_params <- rep(0, ncol(X))
  }
  
  # Estimation des paramètres par optimisation
  opt_result <- optim(par = start_params,
                      fn = log_likelihood,
                      X = X,
                      y = y,
                      method = "BFGS",
                      hessian = TRUE)
  
  estimates <- opt_result$par
  hessian <- opt_result$hessian
  std_errors <- sqrt(diag(solve(hessian)))
  t_values <- estimates / std_errors
  p_values <- 2 * pt(-abs(t_values), df = length(y) - length(estimates))
  
  result <- list(
    coefficients = estimates,
    std_errors = std_errors,
    p_values = p_values,
    log_likelihood = -opt_result$value,
    convergence = opt_result$convergence
  )
  
  return(result)
}

