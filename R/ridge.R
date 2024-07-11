#' Ridge Regression (L2 Regularization)
#'
#' @param X Une matrice de prédicteurs.
#' @param y Un vecteur de la variable dépendante.
#' @param lambda Le paramètre de pénalité.
#' @param tol La tolérance pour la convergence.
#' @param max_iter Le nombre maximal d'itérations.
#' @param verbose Un booléen pour imprimer les informations d'itération.
#' @return Une liste contenant les coefficients estimés, les erreurs standard, les valeurs t, les p-values et les statistiques globales.
#' @export
ridge_regression <- function(X, y, lambda = 0.1, tol = 1e-6, max_iter = 1000, verbose = FALSE) {
  # Fonction pour calculer la perte Ridge
  ridge_loss <- function(w, X, y, lambda) {
    residuals <- y - X %*% w
    loss <- sum(residuals^2) + lambda * sum(w^2)
    return(loss)
  }
  
  # Initialiser les coefficients
  beta_init <- rep(0, ncol(X))
  
  # Optimiser les coefficients en utilisant la méthode BFGS
  fit <- optim(
    par = beta_init,
    fn = ridge_loss,
    X = X,
    y = y,
    lambda = lambda,
    method = 'BFGS',
    control = list(reltol = tol, maxit = max_iter)
  )
  
  beta_ridge <- fit$par
  
  # Calculer les résidus
  residuals <- y - X %*% beta_ridge
  
  # Calculer les erreurs standard des coefficients
  XtX <- t(X) %*% X + lambda * diag(ncol(X))
  XtX_inv <- solve(XtX)
  sigma2 <- sum(residuals^2) / (nrow(X) - ncol(X))
  se_beta_ridge <- sqrt(diag(XtX_inv) * sigma2)
  
  # Calculer les valeurs t et les p-values
  t_values <- beta_ridge / se_beta_ridge
  p_values <- 2 * pt(-abs(t_values), df = nrow(X) - ncol(X))
  
  # Calculer le R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  
  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(X)),
    Estimate = beta_ridge,
    `Std. Error` = se_beta_ridge,
    `t value` = t_values,
    `Pr(>|t|)` = p_values
  )
  
  # Créer un tableau pour les statistiques globales
  global_stats <- data.frame(
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared,
    `Residual Standard Error` = sqrt(sigma2),
    `F-statistic` = NA, # Not applicable for Ridge regression directly
    `F-test P-value` = NA  # Not applicable for Ridge regression directly
  )
  
  # Créer un tableau formaté
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )
  
  return(results)
}
