#' Generalized Least Squares (GLS) Estimation
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#' @param max_iter Maximum number of iterations for the GLS estimation process.
#' @param tol Tolerance level for convergence.
#'
#' @return A list containing the estimated coefficients, standard errors, p-values, adjusted R-squared, F-statistic, and p-value of the F-test.
#' @export
gls_model <- function(data, max_iter = 100, tol = 1e-6) {
  # Vérifier si le data frame contient au moins deux colonnes
  if (ncol(data) < 2) {
    stop("The data frame must contain at least two columns: one dependent variable and at least one independent variable.")
  }
  
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Estimation initiale avec OLS
  beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  residuals <- y - X %*% beta_ols
  
  # Initialiser Sigma comme une matrice identité
  Sigma <- diag(nrow(data))
  
  # Itérer pour ajuster Sigma et les coefficients GLS
  for (i in 1:max_iter) {
    # Estimer la nouvelle Sigma
    Sigma_new <- diag(as.vector(residuals^2))
    
    # Calculer les nouveaux coefficients GLS
    Xt_Sigma_inv <- t(X) %*% solve(Sigma_new)
    beta_gls <- solve(Xt_Sigma_inv %*% X) %*% Xt_Sigma_inv %*% y
    
    # Calculer les nouveaux résidus
    residuals_new <- y - X %*% beta_gls
    
    # Vérifier la convergence
    if (max(abs(residuals_new - residuals)) < tol) {
      residuals <- residuals_new
      Sigma <- Sigma_new
      break
    }
    
    # Mettre à jour les résidus et Sigma
    residuals <- residuals_new
    Sigma <- Sigma_new
  }
  
  # Calculer les erreurs standard des coefficients GLS
  Xt_Sigma_inv_X_inv <- solve(Xt_Sigma_inv %*% X)
  se_beta_gls <- sqrt(diag(Xt_Sigma_inv_X_inv))
  
  # Calculer les statistiques de test robustes
  t_values_robust <- beta_gls / se_beta_gls
  p_values_robust <- 2 * pt(-abs(t_values_robust), df = nrow(X) - ncol(X))
  
  # Calculer le R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  
  # Calculer le F-statistic
  MS_reg <- (SS_tot - SS_res) / (ncol(X) - 1)
  MS_res <- SS_res / (nrow(X) - ncol(X))
  f_statistic <- MS_reg / MS_res
  f_pvalue <- pf(f_statistic, ncol(X) - 1, nrow(X) - ncol(X), lower.tail = FALSE)
  
  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(data)[-1]),
    Value = beta_gls,
    `Standard Error` = se_beta_gls,
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



