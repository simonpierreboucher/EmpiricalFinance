#' White Test for Heteroscedasticity
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#' @param model_results A list containing the results from the `linear_model` function.
#'
#' @return A list containing the White test statistic and the corresponding p-value.
#' @export
white_test <- function(data, model_results) {
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Extraire les résidus du modèle linéaire
  residuals <- y - X %*% as.vector(model_results$coefficients$Value)
  
  # Calculer les résidus au carré
  residuals_squared <- residuals^2
  
  # Construire la matrice des variables indépendantes pour le test de White
  X_white <- cbind(X, X[, -1]^2, apply(X[, -1], 1, prod))
  
  # Ajuster un modèle linéaire pour les résidus au carré en fonction des nouvelles variables indépendantes
  auxiliary_model <- lm(residuals_squared ~ X_white - 1)  # Le "- 1" omet l'intercept puisque X_white a déjà une colonne de 1
  auxiliary_r_squared <- summary(auxiliary_model)$r.squared
  
  # Calculer le test statistique de White
  n <- nrow(data)
  white_statistic <- n * auxiliary_r_squared
  
  # Calculer la p-value associée
  df <- ncol(X_white) - 1  # Degrés de liberté
  white_pvalue <- 1 - pchisq(white_statistic, df = df)
  
  # Créer un tableau formaté
  results <- list(
    white_statistic = white_statistic,
    white_pvalue = white_pvalue
  )
  
  return(results)
}
