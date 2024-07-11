#' Fixed Effects Model Estimation
#'
#' @param data A data frame containing the panel data.
#' @param id_var A string representing the name of the column containing the individual (panel) identifier.
#' @param time_var A string representing the name of the column containing the time variable.
#' @param dep_var A string representing the name of the dependent variable.
#' @param indep_vars A vector of strings representing the names of the independent variables.
#'
#' @return A list containing the estimated coefficients, standard errors, t-values, p-values, and residuals.
#' @export
fixed_effects_model <- function(data, id_var, time_var, dep_var, indep_vars) {
  # Créer les matrices de données
  y <- data[[dep_var]]
  X <- model.matrix(~ . - 1, data = data[, indep_vars])
  
  # Obtenir les identifiants et les temps
  id <- data[[id_var]]
  time <- data[[time_var]]
  
  # Démeaner les données
  y_demeaned <- y - ave(y, id, FUN = mean)
  X_demeaned <- apply(X, 2, function(x) x - ave(x, id, FUN = mean))
  
  # Estimer le modèle par OLS
  beta_fe <- solve(t(X_demeaned) %*% X_demeaned) %*% t(X_demeaned) %*% y_demeaned
  residuals <- y_demeaned - X_demeaned %*% beta_fe
  
  # Calculer les erreurs standard
  N <- length(y)
  K <- ncol(X_demeaned)
  sigma_squared <- sum(residuals^2) / (N - K)
  var_beta_fe <- sigma_squared * solve(t(X_demeaned) %*% X_demeaned)
  se_beta_fe <- sqrt(diag(var_beta_fe))
  
  # Calculer les p-values
  t_values <- beta_fe / se_beta_fe
  p_values <- 2 * pt(-abs(t_values), df = N - K)
  
  # Retourner les résultats
  results <- list(
    coefficients = data.frame(
      Coefficient = colnames(X_demeaned),
      Estimate = beta_fe,
      `Std. Error` = se_beta_fe,
      `t value` = t_values,
      `Pr(>|t|)` = p_values
    ),
    residuals = residuals
  )
  
  return(results)
}
