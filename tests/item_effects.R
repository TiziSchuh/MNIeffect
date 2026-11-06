source("R/raw_functions.R")

nu <- 0
lambda <- .72
nu_g <- 0
lambda_g <- .84
alpha_g <- 0
psi_g <- 1

EDitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)
EADitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)
RESDitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)



# test: compare to numerically approximated integral

D <- function(theta) {
  out <- ((nu_g + lambda_g*theta) - (nu + lambda*theta)) * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}
AD <- function(theta) {
  out <- abs((nu_g + lambda_g*theta) - (nu + lambda*theta)) * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}
SD <- function(theta) {
  out <- ((nu_g + lambda_g*theta) - (nu + lambda*theta))^2 * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}

summary <- data.frame(
  analytical = c(
    EDitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g),
    EADitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g),
    RESDitem(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)),
  numerical = c(
    integrate(D, lower = -10, upper = 10)$value,
    integrate(AD, lower = -10, upper = 10)$value,
    sqrt(integrate(SD, lower = -10, upper = 10)$value) 
  ),
  stringsAsFactors = FALSE
)

print(summary)