#' Linear Model Estimation
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#'
#' @return A list containing the estimated coefficients, standard errors, p-values, adjusted R-squared, F-statistic, and p-value of the F-test.
#' @import stats
#' @export
linear_model <- function(data) {
  # Vérifier si le data frame contient au moins deux colonnes
  if (ncol(data) < 2) {
    stop("The data frame must contain at least two columns: one dependent variable and at least one independent variable.")
  }
  
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Calculer les coefficients avec la pseudo-inverse de Moore-Penrose
  XtX <- t(X) %*% X
  XtX_inv <- tryCatch({
    solve(XtX)
  }, error = function(e) {
    cat("System is computationally singular, using pseudo-inverse.\n")
    ginv(XtX)
  })
  XtY <- t(X) %*% y
  beta <- XtX_inv %*% XtY
  
  # Calculer les valeurs prédictes et les résidus
  y_hat <- X %*% beta
  residuals <- y - y_hat
  
  # Calculer les erreurs standard des coefficients
  n <- nrow(data)
  p <- ncol(X)
  sigma2 <- sum(residuals^2) / (n - p)
  se_beta <- sqrt(diag(XtX_inv) * sigma2)
  
  # Calculer les statistiques de test
  t_values <- beta / se_beta
  p_values <- 2 * pt(-abs(t_values), df = n - p)
  
  # Calculer le R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p))
  
  # Calculer le F-statistic et la p-value du F-test
  MS_reg <- (SS_tot - SS_res) / (p - 1)
  MS_res <- SS_res / (n - p)
  f_statistic <- MS_reg / MS_res
  f_pvalue <- pf(f_statistic, p - 1, n - p, lower.tail = FALSE)
  
  # Créer un tableau formaté pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(data)[-1]),
    Value = beta,
    `Standard Error` = se_beta,
    `P-value` = p_values
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

#' Generalized Inverse of a Matrix (Moore-Penrose)
#'
#' @param A A matrix to be inverted.
#'
#' @return The Moore-Penrose pseudo-inverse of matrix A.
ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  d <- s$d
  d[d > tol] <- 1 / d[d > tol]
  d[d <= tol] <- 0
  s$v %*% diag(d) %*% t(s$u)
}


