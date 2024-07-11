#' LASSO (Least Absolute Shrinkage and Selection Operator) Regression
#'
#' @param X Une matrice de prédicteurs.
#' @param y Un vecteur de la variable dépendante.
#' @param lambda Le paramètre de pénalité.
#' @param tol La tolérance pour la convergence.
#' @param max_iter Le nombre maximal d'itérations.
#' @param verbose Un booléen pour imprimer les informations d'itération.
#' @return Un vecteur contenant les coefficients estimés.
#' @export
lasso_regression <- function(X, y, lambda = 0.1, tol = 1e-6, max_iter = 1000, verbose = FALSE) {
  # Fonction de seuil soft
  soft_thresh <- function(a, b) {
    sign(a) * pmax(0, abs(a) - b)
  }
  
  # Initialiser les coefficients
  beta <- solve(t(X) %*% X + diag(lambda, ncol(X))) %*% t(X) %*% y
  tol_curr <- tol + 1
  iter <- 1
  
  while (tol_curr > tol && iter <= max_iter) {
    beta_old <- beta
    
    for (j in 1:ncol(X)) {
      # Calculer la somme des autres coefficients
      X_j <- X[, j]
      r_j <- y - X %*% beta + X_j * beta[j]
      beta[j] <- soft_thresh(sum(X_j * r_j), lambda) / sum(X_j^2)
    }
    
    tol_curr <- sum((beta - beta_old)^2)
    
    if (verbose && iter %% 10 == 0) {
      cat("Iteration:", iter, "Tolerance:", tol_curr, "\n")
    }
    
    iter <- iter + 1
  }
  
  return(beta)
}
