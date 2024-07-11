#' Maximum Likelihood Estimation (MLE) for Linear Model
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#'
#' @return A list containing the estimated coefficients, standard errors, p-values, adjusted R-squared, log-likelihood, AIC, and BIC.
#' @export
mle_model <- function(data) {
  # Vérifier si le data frame contient au moins deux colonnes
  if (ncol(data) < 2) {
    stop("The data frame must contain at least two columns: one dependent variable and at least one independent variable.")
  }
  
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Définir la fonction de log-vraisemblance
  log_likelihood <- function(params) {
    beta <- params[1:ncol(X)]
    sigma <- exp(params[length(params)])  # Assurer que sigma est positif
    residuals <- y - X %*% beta
    ll <- -0.5 * nrow(X) * log(2 * pi * sigma^2) - (1 / (2 * sigma^2)) * sum(residuals^2)
    return(-ll)  # On minimise donc on retourne le négatif
  }
  
  # Initialiser les paramètres (beta et log(sigma))
  init_params <- c(rep(0, ncol(X)), log(sd(y)))
  
  # Optimisation par maximum de vraisemblance
  mle_optim <- optim(par = init_params, fn = log_likelihood, method = "BFGS", hessian = TRUE)
  
  # Récupérer les paramètres optimisés
  beta_mle <- mle_optim$par[1:ncol(X)]
  sigma_mle <- exp(mle_optim$par[length(mle_optim$par)])
  log_lik <- -mle_optim$value
  
  # Calculer les erreurs standard des coefficients MLE
  fisher_info <- solve(mle_optim$hessian)
  se_params <- sqrt(diag(fisher_info))
  se_beta_mle <- se_params[1:ncol(X)]
  
  # Calculer les statistiques de test
  t_values <- beta_mle / se_beta_mle
  p_values <- 2 * pt(-abs(t_values), df = nrow(X) - ncol(X))
  
  # Calculer le R-carré et le R-carré ajusté
  SS_tot <- sum((y - mean(y))^2)
  residuals <- y - X %*% beta_mle
  SS_res <- sum(residuals^2)
  r_squared <- 1 - (SS_res / SS_tot)
  adj_r_squared <- 1 - (1 - r_squared) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  
  # Calculer AIC et BIC
  aic <- 2 * length(init_params) - 2 * log_lik
  bic <- log(nrow(X)) * length(init_params) - 2 * log_lik
  
  # Créer un tableau pour les coefficients
  coefficients_table <- data.frame(
    Coefficient = c("(Intercept)", colnames(data)[-1]),
    Value = beta_mle,
    `Standard Error` = se_beta_mle,
    `P-value` = p_values
  )
  rownames(coefficients_table) <- NULL
  
  # Créer un tableau pour les statistiques globales
  global_stats <- data.frame(
    `R-squared` = r_squared,
    `Adjusted R-squared` = adj_r_squared,
    `Log-likelihood` = log_lik,
    `AIC` = aic,
    `BIC` = bic
  )
  
  # Créer un tableau formaté
  results <- list(
    coefficients = coefficients_table,
    global_stats = global_stats
  )
  
  return(results)
}
