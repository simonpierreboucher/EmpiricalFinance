#' Breusch-Pagan Test for Heteroscedasticity
#'
#' @param data A data frame where the first column is the dependent variable (Y) and the remaining columns are the independent variables (X1, X2, ...).
#' @param model_results A list containing the results from the `linear_model` function.
#'
#' @return A list containing the Breusch-Pagan test statistic and the corresponding p-value.
#' @export
breusch_pagan_test <- function(data, model_results) {
  # Extraire les variables dépendantes et indépendantes
  y <- as.matrix(data[, 1])
  X <- as.matrix(cbind(1, data[, -1]))  # Ajouter une colonne de 1 pour l'intercept
  
  # Extraire les résidus du modèle linéaire
  residuals <- y - X %*% as.vector(model_results$coefficients$Value)
  
  # Calculer les résidus au carré
  residuals_squared <- residuals^2
  
  # Ajuster un modèle linéaire pour les résidus au carré en fonction des variables indépendantes
  auxiliary_model <- lm(residuals_squared ~ X - 1)  # Le "- 1" omet l'intercept puisque X a déjà une colonne de 1
  auxiliary_r_squared <- summary(auxiliary_model)$r.squared
  
  # Calculer le test statistique de Breusch-Pagan
  n <- nrow(data)
  bp_statistic <- n * auxiliary_r_squared
  
  # Calculer la p-value associée
  df <- ncol(X) - 1  # Degrés de liberté
  bp_pvalue <- 1 - pchisq(bp_statistic, df = df)
  
  # Créer un tableau formaté
  results <- list(
    bp_statistic = bp_statistic,
    bp_pvalue = bp_pvalue
  )
  
  return(results)
}
