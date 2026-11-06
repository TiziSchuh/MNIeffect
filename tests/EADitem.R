source("R/raw_functions.R")

nu <- 0
lambda <- .72
nu_g <- 0
lambda_g <- .84
alpha_g <- 0
psi_g <- 1

EADitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)


# test: compare to numerically approximated integral

dRf <- function(theta) {
  out <- abs((nu_g + lambda_g*theta) - (nu + lambda*theta)) * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}

integrate(dRf, lower = -10, upper = 10) 