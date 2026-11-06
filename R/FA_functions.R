#
# implementation based on:
# 
# Schuhbeck, T. M. B., Sterner, P., & Goretzko, D. (2026). 
# Quantifying Measurement Non-Invariance Beyond Simple Structure: 
# The Closed Formulas of Universal Effect Size Measures for MI. 
# Structural Equation Modeling: A Multidisciplinary Journal, 33(1), 81-88.
#
# =============================================================================


library(pracma) # erf() needed in \code(exp_response_diff)


exp_response_diff <- function(nu, lambda, nu_g, lambda_g, alpha_g, psi_g,
                              difference = "absolute") {
  # nu: intercept
  # lambda: loading(-vector)
  # alpha: latent mean(-vector)
  # psi: latent variance(/covariance matrix)
  
  dnu <- nu_g - nu
  dlambda <- matrix(lambda_g) - matrix(lambda)

  # dR: expected response difference
  dR_mean <- dnu + t(dlambda) %*% matrix(alpha_g)
  dR_var <- t(dlambda) %*% psi_g %*% dlambda

  if (difference == "signed") {
    exp_dR <- dR_mean
  } else if (difference == "absolute") {
    exp_dR <- (sqrt(dR_var) * sqrt(2/pi) * exp(-((dR_mean^2) / (2*dR_var)))
            + dR_mean * erf(dR_mean / sqrt(2 * dR_var)))
    if (dR_mean == 0 & dR_var == 0) {exp_dR <- matrix(0)}
  } else if (difference == "squared") {
    exp_dR <- dR_mean^2 + dR_var
  } else {
    stop("only \"signed\", \"absolute\" and \"squared\" differences supported")
  }

  exp_dR[1, 1] # get rid of matrix
}


EDitem <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(nu, lamda, nu_g, lambda_g, alpha_g, psi_g, 
                    difference = "signed")
}

EADitem <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(nu, lamda, nu_g, lambda_g, alpha_g, psi_g, 
                    difference = "absolute")
}

ESDitem <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(nu, lamda, nu_g, lambda_g, alpha_g, psi_g, 
                    difference = "squared")
}

RESDitem <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  sqrt(ESDitem(nu, lamda, nu_g, lambda_g, alpha_g, psi_g))
} 

EDtest <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(sum(nu), sum(lamda), sum(nu_g), sum(lambda_g), alpha_g, psi_g, 
                    difference = "signed")
}

EADtest <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(sum(nu), sum(lamda), sum(nu_g), sum(lambda_g), alpha_g, psi_g, 
                    difference = "absolute")
}

ESDtest <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  exp_response_diff(sum(nu), sum(lamda), sum(nu_g), sum(lambda_g), alpha_g, psi_g, 
                    difference = "squared")
}

RESDtest <- function(nu, lamda, nu_g, lambda_g, alpha_g, psi_g) {
  sqrt(ESDtest(nu, lamda, nu_g, lambda_g, alpha_g, psi_g))
} 
