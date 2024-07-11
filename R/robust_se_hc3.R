#' Calculate Robust Standard Errors (HC3)
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#' @param model_results A list containing the results from the `linear_model` function.
#'
#' @return A list containing the estimated coefficients, robust standard errors, p-values, adjusted R-squared, F-statistic, and p-value of the F-test.
#' @export
robust_se_hc3 <- function(data, model_results) {
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Extraire les coefficients estimés
  beta <- as.vector(model_results$coefficients$Value)
  
  # Calculer les valeurs prédictes et les résidus
  y_hat <- X %*% beta
  residuals <- y - y_hat
  
  # Calculer la matrice de covariances robuste HC3
  n <- nrow(data)
  p <- ncol(X)
  X_t <- t(X)
  leverage <- diag(X %*% solve(X_t %*% X) %*% X_t)
  residuals_adjusted <- residuals^2 / (1 - leverage)^2
  residuals_diag <- diag(as.vector(residuals_adjusted))
  meat <- X_t %*% residuals_diag %*% X
  sandwich <- solve(X_t %*% X) %*% meat %*% solve(X_t %*% X)
  
  # Erreurs standard robustes
  diag_sandwich <- diag(sandwich)
  diag_sandwich[diag_sandwich < 0] <- 0  # Ajuster les valeurs négatives
  robust_se <- sqrt(diag_sandwich)
  
  # Calculer les statistiques de test robustes
  t_values_robust <- beta / robust_se
  p_values_robust <- 2 * pt(-abs(t_values_robust), df = n - p)
  
  # Calculer le R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p))
  
  # Calculer le F-statistic
  MS_reg <- (SS_tot - SS_res) / (p - 1)
  MS_res <- SS_res / (n - p)
  f_statistic <- MS_reg / MS_res
  f_pvalue <- pf(f_statistic, p - 1, n - p, lower.tail = FALSE)
  
  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(data)[-1]),
    Value = beta,
    `Robust Standard Error` = robust_se,
    `P-value` = p_values_robust
  )
  rownames(coefficients_table) <- NULL
  
  # Créer un tableau pour les statistiques globales
  global_stats <- data.frame(
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared,
    `F-statistic` = f_statistic,
    `F-test P-value` = f_pvalue
  )
  
  # Créer un tableau formaté
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )
  
  return(results)
}

