source("R/FA_functions.R")


MNIeffect <- function(fit, baseline, std = NULL, sumscore.include = TRUE, measure = "EAD") {
  
  lambda <- lapply(lavInspect(fit, "est"), function(x) {unname(x$lambda)})
  nu <- lapply(lavInspect(fit, "est"), function(x) {unname(x$nu)})
  alpha <- lapply(lavInspect(fit, "est"), function(x) {unname(x$alpha)})
  psi <- lapply(lavInspect(fit, "est"), function(x) {unname(x$psi)})

  nitems <- length(nu[[1]])
  ngroups <- lavInspect(fit, "ngroups")
  nobs <- lavInspect(fit, "nobs")
  ntotal <- lavInspect(fit, "ntotal")
  p <- nobs/ntotal # group proportions

  if (class(baseline) == "lavaan") {
    lambda_base <- lavTech(baseline, "est")$lambda
    
    if (is.null(lavTech(baseline, "est")$nu)){
      stop("'meanstructure = TRUE' required in the baseline fit object")
    }
    nu_base <- lavTech(baseline, "est")$nu
  }
  else if (baseline == "pooled") {
    lambda_base <- 0
    nu_base <- 0
    for (g in 1:ngroups) {
      lambda_base <- lambda_base + p[g] * lambda[[g]]
      nu_base <- nu_base + p[g] * nu[[g]]
    }
  }
  else {
    lambda_base <- unname(lambda[baseline])[[1]]
    nu_base <- unname(nu[baseline])[[1]]
  }

  # ToDo: standardization
  # cov <- lavTech(fitGlobal, "sampstat")[[1]]$cov
  # sd <- sqrt(diag(cov))

  
  effect <- matrix(0, ngroups, nitems)
  effect_total <- numeric(nitems)
  effect_sumscore <- matrix(0, ngroups, ncol = 1)
  effect_sumscore_total <- 0

  if(measure == "ED"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- EDitem(nu_base[i], lambda_base[i, ], 
                              nu[[g]][i], lambda[[g]][i, ], 
                              alpha[[g]], psi[[g]])
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- EDtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]])
      effect_sumscore_total <- effect_sumscore_total + p[g] * effect_sumscore[g]
    }
  }

  if(measure == "EAD"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- EADitem(nu_base[i], lambda_base[i, ], 
                               nu[[g]][i], lambda[[g]][i, ], 
                               alpha[[g]], psi[[g]])
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- EADtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]])
      effect_sumscore_total <- effect_sumscore_total + p[g] * effect_sumscore[g]
    }
  }

  if(measure == "ESD"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- ESDitem(nu_base[i], lambda_base[i, ], 
                              nu[[g]][i], lambda[[g]][i, ], 
                              alpha[[g]], psi[[g]])
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- ESDtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]])
      effect_sumscore_total <- effect_sumscore_total + p[g] * effect_sumscore[g]
    }
  }

  if(measure == "RESD"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- ESDitem(nu_base[i], lambda_base[i, ], 
                              nu[[g]][i], lambda[[g]][i, ], 
                              alpha[[g]], psi[[g]])
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
        effect[g,i] <- sqrt(effect[g,i])
      }
      effect_total[i] <- sqrt(effect_total[i])
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- ESDtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]])
      effect_sumscore_total <- effect_sumscore_total + p[g] * effect_sumscore[g]
      effect_sumscore[g] <- sqrt(effect_sumscore[g])
    }
    effect_sumscore_total <- sqrt(effect_sumscore_total)
  }

  frame <- data.frame(
    rbind(effect, effect_total),
    row.names = append(lavInspect(fit, what = "group.label"), "Total")
  )
  if (sumscore.include == TRUE) {
    frame_sumscore <- data.frame(
      `Sum-Score` = rbind(effect_sumscore, effect_sumscore_total)
    )
    frame <- cbind(frame,frame_sumscore)
  }

  frame
}



#' Expected Difference
#' 
#' This function calculates the Expected Difference (ED) and Total Expected
#' Difference (TED) of responses for linear common factor models with 
#' (jointly) normally distributed factors  
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of th reference group as a string.
#' @param std Default `NULL`.   
#' 
#' @examples
#' ED(fit, baseline = "pooled")
#' ED(fit, baseline = fit_global)
#' ED(fit, baseline = "women")
#' 
#' @returns All item-level values ED (for each group) and TED (across groups) 
#' @export
ED <- function(fit, baseline = "pooled", std = NULL) {
  MNIeffect(fit, baseline, std, measure = "ED")
}

EAD <- function(fit, baseline = "pooled", std = NULL) {
  MNIeffect(fit, baseline, std, measure = "EAD")
}

ESD <- function(fit, baseline = "pooled", std = NULL) {
  MNIeffect(fit, baseline, std, measure = "ESD")
}

RESD <- function(fit, baseline = "pooled", std = NULL) {
  MNIeffect(fit, baseline, std, measure = "RESD")
}