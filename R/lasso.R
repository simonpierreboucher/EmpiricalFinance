#' LASSO Regression (L1 Regularization)
#'
#' @param X Une matrice de prédicteurs.
#' @param y Un vecteur de la variable dépendante.
#' @param lambda Le paramètre de pénalité.
#' @param tol La tolérance pour la convergence.
#' @param max_iter Le nombre maximal d'itérations.
#' @param verbose Un booléen pour imprimer les informations d'itération.
#' @return Une liste contenant les coefficients estimés, les erreurs standard, les valeurs t, les p-values et les statistiques globales.
#' @export
lasso_regression <- function(X, y, lambda = 0.1, tol = 1e-6, max_iter = 1000, verbose = FALSE) {
  # Fonction pour le seuillage doux
  soft_thresh <- function(a, b) {
    out <- rep(0, length(a))
    out[a > b] <- a[a > b] - b
    out[a < -b] <- a[a < -b] + b
    out
  }
  
  # Initialiser les coefficients
  beta <- rep(0, ncol(X))
  tol_curr <- 1
  iter <- 1
  N <- nrow(X)
  J <- ncol(X)
  
  # Itération pour les coefficients LASSO
  while (tol_curr > tol && iter < max_iter) {
    beta_old <- beta
    for (j in 1:J) {
      rho <- sum(X[, j] * (y - X %*% beta + X[, j] * beta[j]))
      beta[j] <- soft_thresh(rho, lambda * N) / sum(X[, j]^2)
    }
    tol_curr <- sum((beta - beta_old)^2)
    iter <- iter + 1
    if (verbose && iter %% 10 == 0) message("Iteration: ", iter, " Tolerance: ", tol_curr)
  }
  
  # Calcul des résidus
  residuals <- y - X %*% beta
  
  # Calcul des erreurs standard approximatives par permutation
  n_permutations <- 100
  beta_permuted <- matrix(0, n_permutations, J)
  for (i in 1:n_permutations) {
    y_permuted <- sample(y)
    beta_perm <- rep(0, J)
    tol_curr_perm <- 1
    iter_perm <- 1
    while (tol_curr_perm > tol && iter_perm < max_iter) {
      beta_old_perm <- beta_perm
      for (j in 1:J) {
        rho_perm <- sum(X[, j] * (y_permuted - X %*% beta_perm + X[, j] * beta_perm[j]))
        beta_perm[j] <- soft_thresh(rho_perm, lambda * N) / sum(X[, j]^2)
      }
      tol_curr_perm <- sum((beta_perm - beta_old_perm)^2)
      iter_perm <- iter_perm + 1
    }
    beta_permuted[i, ] <- beta_perm
  }
  
  se_beta <- apply(beta_permuted, 2, sd)
  
  # Calcul des valeurs t et des p-values
  t_values <- beta / se_beta
  p_values <- 2 * pt(-abs(t_values), df = N - J)
  
  # Calcul du R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((N - 1) / (N - J))
  
  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(X)),
    Estimate = beta,
    `Std. Error` = se_beta,
    `t value` = t_values,
    `Pr(>|t|)` = p_values
  )
  
  # Créer un tableau pour les statistiques globales
  global_stats <- data.frame(
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared,
    `Residual Standard Error` = sqrt(sum(residuals^2) / (N - J)),
    `F-statistic` = NA, # Not applicable for LASSO directly
    `F-test P-value` = NA  # Not applicable for LASSO directly
  )
  
  # Créer un tableau formaté
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )
  
  return(results)
}

