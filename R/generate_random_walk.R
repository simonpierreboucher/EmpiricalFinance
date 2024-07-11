#' Générer une série temporelle de type Random Walk
#'
#' @param n Un entier représentant le nombre de périodes.
#' @param drift Un nombre représentant la dérive (drift). Par défaut à 0.
#' @param trend Un vecteur de tendances. Par défaut à NULL (pas de tendance).
#'
#' @return Un vecteur numérique représentant la série temporelle générée.
#' @export
generate_random_walk <- function(n, drift = 0, trend = NULL) {
  # Générer des innovations aléatoires
  innovations <- rnorm(n)
  
  # Initialiser la série
  series <- numeric(n)
  series[1] <- innovations[1]
  
  # Générer la série
  for (i in 2:n) {
    if (!is.null(trend)) {
      series[i] <- series[i - 1] + innovations[i] + drift + trend[i]
    } else {
      series[i] <- series[i - 1] + innovations[i] + drift
    }
  }
  
  return(series)
}

