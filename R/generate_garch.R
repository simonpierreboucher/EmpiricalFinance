#' Générer une série temporelle de type GARCH(p, q)
#'
#' @param n Un entier représentant le nombre de périodes.
#' @param omega Un nombre représentant le terme constant du modèle GARCH.
#' @param alpha Un vecteur de coefficients ARCH.
#' @param beta Un vecteur de coefficients GARCH.
#'
#' @return Un vecteur numérique représentant la série temporelle générée.
#' @export
generate_garch <- function(n, omega, alpha, beta) {
  p <- length(alpha)
  q <- length(beta)
  
  # Initialiser les innovations et les variances conditionnelles
  innovations <- numeric(n + max(p, q))
  conditional_variances <- numeric(n + max(p, q))
  
  # Initialiser la série
  series <- numeric(n + max(p, q))
  
  # Générer la série
  for (i in (max(p, q) + 1):(n + max(p, q))) {
    conditional_variance <- omega + sum(alpha * innovations[(i - p):(i - 1)]^2) + sum(beta * conditional_variances[(i - q):(i - 1)])
    conditional_variances[i] <- conditional_variance
    innovations[i] <- rnorm(1, mean = 0, sd = sqrt(conditional_variance))
    series[i] <- innovations[i]
  }
  
  return(series[(max(p, q) + 1):(n + max(p, q))])
}

# Exemple d'utilisation
garch_series <- generate_garch(1000, omega = 0.1, alpha = c(0.1, 0.2), beta = c(0.3))
plot(garch_series, type = "l", main = "GARCH(2, 1) Series", xlab = "Time", ylab = "Value")
