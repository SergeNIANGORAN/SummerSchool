
#' Title : Creation of the normal multivariate distribution function
#' Description of the function
#' This function compute for a univariable or mutivariable parameter
#'  the PDF (Probability Distribution Function) for a normal density
#'
#'  The details
#'
#'
#'
#' @param x a p x n data matrix with n the number of observations and
#'p the number of dimensions
#' @param mean mean vector
#' @param varcovM variance-covariance matrix
#' @param Log logical flag for returning the log of the probability density
#'function. Default is \code{TRUE}.
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @return a list containing the input matrix x and y the multivariate-Normal probability density function
#' computed at x
#'
#' @export
#'
#' @examples
#' mvnpdf(x=matrix(1.96), Log=FALSE)
#'dnorm(1.96)
#'
#'mvnpdf(x=matrix(rep(1.96, 2), nrow=2, ncol=1), Log=FALSE)
#'
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  # On peut utiliser differentes options pour l'initialisation
  # y <- vector(mode="numeric", length = n)
  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}


# Fonction PDF
