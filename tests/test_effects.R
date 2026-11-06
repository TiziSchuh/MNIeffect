source("R/FA_functions.R")

# ------------------------------------------------------------------------------------------------------
# three items, one factor

nu <- c(0, 0, 0)
lambda <- c(.5, .6, .7)
nu_g <- c(0, 0, 0)
lambda_g <- c(.6, .7, .8)
alpha_g <- 0
psi_g <- 1

EDtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)
EADtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)
RESDtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)


# ------------------------------------------------------------------------------------------------------
# compare to numerically approximated integral

diff <- function(theta) {
  out <- (sum(nu_g + lambda_g*theta) - (nu + lambda*theta)) * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}
diff_abs <- function(theta) {
  out <- abs(sum((nu_g + lambda_g*theta) - (nu + lambda*theta))) * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}
diff_sqr <- function(theta) {
  out <- sum((nu_g + lambda_g*theta) - (nu + lambda*theta))^2 * dnorm(theta, mean=alpha_g, sd=sqrt(psi_g))
}

summary <- data.frame(
  analytical = c(
    EDtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g),
    EADtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g),
    RESDtest(nu, lambda, nu_g, lambda_g, alpha_g, psi_g)),
  numerical = c(
    integrate(diff, lower = -10, upper = 10)$value,
    integrate(diff_abs, lower = -10, upper = 10)$value,
    sqrt(integrate(diff_sqr, lower = -10, upper = 10)$value) 
  ),
  stringsAsFactors = FALSE
)

print(summary)
