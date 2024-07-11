#' Random Effects Model Estimation
#'
#' @param data A data frame containing the panel data.
#' @param id_var A string representing the name of the column containing the individual (panel) identifier.
#' @param time_var A string representing the name of the column containing the time variable.
#' @param dep_var A string representing the name of the dependent variable.
#' @param indep_vars A vector of strings representing the names of the independent variables.
#'
#' @return A list containing the estimated coefficients, standard errors, t-values, p-values, and residuals.
#' @export
random_effects_model <- function(data, id_var, time_var, dep_var, indep_vars) {
  # Créer les matrices de données
  y <- data[[dep_var]]
  X <- model.matrix(~ . - 1, data = data[, indep_vars])
  
  # Obtenir les identifiants et les temps
  id <- data[[id_var]]
  time <- data[[time_var]]
  
  # Estimer le modèle par OLS pour obtenir les résidus
  ols_model <- lm(as.formula(paste(dep_var, "~", paste(indep_vars, collapse = " + "))), data = data)
  ols_residuals <- residuals(ols_model)
  
  # Calculer la variance des résidus et la variance des effets aléatoires
  n <- length(unique(id))
  T <- length(unique(time))
  sigma_u2 <- sum((ave(ols_residuals, id, FUN = mean))^2) / (n - 1)
  sigma_e2 <- sum((ols_residuals - ave(ols_residuals, id, FUN = mean))^2) / (n * (T - 1))
  theta <- 1 - sqrt(sigma_e2 / (sigma_e2 + T * sigma_u2))
  
  # Démeaner les données en utilisant theta
  y_demeaned <- y - theta * ave(y, id, FUN = mean)
  X_demeaned <- apply(X, 2, function(x) x - theta * ave(x, id, FUN = mean))
  
  # Estimer le modèle par OLS sur les données démeanées
  beta_re <- solve(t(X_demeaned) %*% X_demeaned) %*% t(X_demeaned) %*% y_demeaned
  residuals <- y_demeaned - X_demeaned %*% beta_re
  
  # Calculer les erreurs standard
  N <- length(y)
  K <- ncol(X_demeaned)
  sigma_squared <- sum(residuals^2) / (N - K)
  var_beta_re <- sigma_squared * solve(t(X_demeaned) %*% X_demeaned)
  se_beta_re <- sqrt(diag(var_beta_re))
  
  # Calculer les p-values
  t_values <- beta_re / se_beta_re
  p_values <- 2 * pt(-abs(t_values), df = N - K)
  
  # Retourner les résultats
  results <- list(
    coefficients = data.frame(
      Coefficient = colnames(X_demeaned),
      Estimate = beta_re,
      `Std. Error` = se_beta_re,
      `t value` = t_values,
      `Pr(>|t|)` = p_values
    ),
    residuals = residuals
  )
  
  return(results)
}
