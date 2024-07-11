#' Générer une série temporelle de type ARIMA
#'
#' @param n Un entier représentant le nombre de périodes.
#' @param ar Un vecteur de coefficients AR.
#' @param ma Un vecteur de coefficients MA.
#' @param d Un entier représentant le nombre de différenciations.
#' @param drift Un nombre représentant la dérive (drift). Par défaut à 0.
#'
#' @return Un vecteur numérique représentant la série temporelle générée.
#' @export
generate_arima <- function(n, ar = NULL, ma = NULL, d = 0, drift = 0) {
  # Générer des innovations aléatoires
  innovations <- rnorm(n + max(length(ar), length(ma), d))
  
  # Initialiser la série
  series <- numeric(n + max(length(ar), length(ma), d))
  
  # Générer la série
  for (i in (max(length(ar), length(ma), d) + 1):(n + max(length(ar), length(ma), d))) {
    ar_part <- if (!is.null(ar)) sum(ar * rev(series[(i - length(ar)):(i - 1)])) else 0
    ma_part <- if (!is.null(ma)) sum(ma * rev(innovations[(i - length(ma)):(i - 1)])) else 0
    series[i] <- ar_part + ma_part + innovations[i] + drift
  }
  
  # Différenciation
  if (d > 0) {
    for (i in 1:d) {
      series <- diff(series)
    }
  }
  
  return(series[(max(length(ar), length(ma), d) + 1):(n + max(length(ar), length(ma), d))])
}

