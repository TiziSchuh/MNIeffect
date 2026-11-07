source("R/FA_functions.R")


MNIeffect <- function(fit, baseline, std = NULL, measure = "EAD") {
  
  lambda <- lapply(lavInspect(fit, "est"), function(x) {
    unname(x$lambda)
  })
  nu <- lapply(lavInspect(fit, "est"), function(x) {
    unname(x$nu)
  })
  alpha <- lapply(lavInspect(fit, "est"), function(x) {
    unname(x$alpha)
  })
  psi <- lapply(lavInspect(fit, "est"), function(x) {
    unname(x$psi)
  })

  if (class(baseline) == "lavaan") {
    lambda_base <- lavTech(baseline, "est")$lambda
    
    if (is.null(lavTech(baseline, "est")$nu)){
      stop("'meanstructure = TRUE' required in the baseline fit object")
    }
    nu_base <- lavTech(baseline, "est")$nu
  }
  else {
    lambda_base <- unname(lambda[baseline])[[1]]
    nu_base <- unname(nu[baseline])[[1]]
  }

  nitems <- length(nu[[1]])
  ngroups <- lavInspect(fit, "ngroups")
  
  nobs <- lavInspect(fit, "nobs")
  ntotal <- lavInspect(fit, "ntotal")
  p <- nobs/ntotal # group proportions

  # ToDo: standardization
  # cov <- lavTech(fitGlobal, "sampstat")[[1]]$cov
  # sd <- sqrt(diag(cov))

  
  effect <- matrix(0, ngroups, nitems)
  effect_total <- numeric(nitems)

  if(measure == "ED"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- EDitem(nu_base[i], lambda_base[i, ], 
                              nu[[g]][i], lambda[[g]][i, ], 
                              alpha[[g]], psi[[g]])
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
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
  }


  frame <- data.frame(
    rbind(effect, effect_total),
    row.names = append(lavInspect(fit, what = "group.label"), "Total")
  )
  
  frame
  }


ED <- function(fit, baseline, std = NULL) {
  MNIeffect(fit, baseline, std, measure = "ED")
}

EAD <- function(fit, baseline, std = NULL) {
  MNIeffect(fit, baseline, std, measure = "EAD")
}

ESD <- function(fit, baseline, std = NULL) {
  MNIeffect(fit, baseline, std, measure = "ESD")
}

RESD <- function(fit, baseline, std = NULL) {
  MNIeffect(fit, baseline, std, measure = "RESD")
}