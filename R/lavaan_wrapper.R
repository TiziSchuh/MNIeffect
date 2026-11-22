source("R/FA_functions.R")


# helper function
get_sumscores <- function(fit) {
  X <- lavInspect(fit, "data")
  if (!is.list(X)) X <- list(X)
  glabs <- lavInspect(fit, "group.label")
  names(X) <- if (length(glabs)) glabs 

  # all observed variables except covariates
  items <- lavNames(fit, "ov.nox")   

  lapply(X, function(M) {
    M <- as.matrix(M)[, items, drop = FALSE]
    rowSums(M[complete.cases(M), , drop = FALSE])        
  })
}


MNIeffect <- function(fit, baseline, std = FALSE, sumscore.include = TRUE, measure = "EAD") {
  
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


  sd <- matrix(1, ngroups, nitems)
  sd_sumscore <- rep(1, ngroups)

  if (std == "global") {
    if (class(baseline) == "lavaan") { 
      sd_base <- sqrt(diag(lavTech(baseline, "sampstat")[[1]]$cov))
      sd <- t(replicate(ngroups, sd_base))
      sd_sumscore <- replicate(ngroups,sd(get_sumscores(baseline)[[1]]))
    }
    else {
      sd_pool <- rep(0, nitems)
      sd_sumscore_pool <- 0 
      for (g in 1:ngroups) {
        sd_pool <- sd_pool + (nobs[g]-1) * diag(lavTech(fit, "sampstat")[[g]]$cov)
        sd_sumscore_pool <- sd_sumscore_pool + (nobs[g]-1) * (sd(get_sumscores(fit)[[g]])^2)
      }
      sd_pool <- sqrt(sd_pool/(sum(nobs)-ngroups))
      sd_sumscore_pool <- sqrt(sd_sumscore_pool/(sum(nobs)-ngroups))
      sd <- t(replicate(ngroups, sd_pool))
      sd_sumscore <- replicate(ngroups, sd_sumscore_pool)
    }
  }
  else if (std == "group") {
    for (g in 1:ngroups) {
      sd[g,] <- sqrt(diag(lavTech(fit, "sampstat")[[g]]$cov))
      sd_sumscore[g] <- sd(get_sumscores(fit)[[g]])
    }
  }
  else if (std == "reference") {
    sd_ref <- sqrt(diag(lavTech(fit, "sampstat")[[baseline]]$cov))
    sd <- t(replicate(ngroups, sd_ref))
    sd_sumscore_ref <- sd(get_sumscores(fit)[[baseline]])
    sd_sumscore <- replicate(ngroups, sd_sumscore_ref)
  }
  else if (std == "pooled") {
    sd_ref <- sqrt(diag(lavTech(fit, "sampstat")[[baseline]]$cov))
    sd_sumscore_ref <- sd(get_sumscores(fit)[[baseline]])
    r <- match(baseline, lavInspect(fit, "group.label"))
    for (g in 1:ngroups) {
      sd[g,] <- (nobs[r]-1)*sd_ref + (nobs[g]-1)*diag(lavTech(fit, "sampstat")[[g]]$cov)
      sd[g,] <- sqrt(sd[g,]/(nobs[g]+nobs[r]-2))
      sd_sumscore[g] <- (nobs[r]-1)*sd_sumscore_ref + (nobs[g]-1)*(sd(get_sumscores(fit)[[g]])^2)
      sd_sumscore[g] <- sqrt(sd_sumscore[g]/(nobs[g]+nobs[r]-2))
    }
  }



  effect <- matrix(0, ngroups, nitems)
  effect_total <- numeric(nitems)
  effect_sumscore <- matrix(0, ngroups, ncol = 1)
  effect_sumscore_total <- 0

  if(measure == "ED"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- EDitem(nu_base[i], lambda_base[i, ], 
                              nu[[g]][i], lambda[[g]][i, ], 
                              alpha[[g]], psi[[g]]) / sd[g,i]
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- EDtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]]) / sd_sumscore[g]
      effect_sumscore_total <- effect_sumscore_total + p[g] * effect_sumscore[g]
    }
  }

  if(measure == "EAD"){ 
    for (i in 1:nitems) {
      for (g in 1:ngroups) {
        effect[g,i] <- EADitem(nu_base[i], lambda_base[i, ], 
                               nu[[g]][i], lambda[[g]][i, ], 
                               alpha[[g]], psi[[g]]) / sd[g,i]
        effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
      }
    }
    for (g in 1:ngroups) {
      effect_sumscore[g] <- EADtest(nu_base, lambda_base, 
                                   nu[[g]], lambda[[g]], 
                                   alpha[[g]], psi[[g]]) / sd_sumscore[g]
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



#' Average Difference of Expected Responses
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
#' @param std Default `FALSE` returns ED in orginial metric. 
#' Option `"global"` returns ED in global/globally pooled SDs. 
#' Option `"group"` returns ED in respective group's SD.
#' Option `"reference"` returns ED in reference group's SD.
#' Option `"pooled"` returns ED in reference-focal-pooled SD.
#' 
#' @examples
#' ED(fit, baseline = "pooled")
#' ED(fit, baseline = fit_global)
#' ED(fit, baseline = "women")
#' 
#' @returns All item-level values ED (for each group) and TED (across groups) 
#' @export
ED <- function(fit, baseline = "pooled", std = FALSE) {
  MNIeffect(fit, baseline, std, measure = "ED")
}

EAD <- function(fit, baseline = "pooled", std = FALSE) {
  MNIeffect(fit, baseline, std, measure = "EAD")
}

MNIdir <- function(fit, baseline = "pooled") {
  df <- ED(fit, baseline, std = FALSE) / EAD(fit, baseline, std = FALSE) 
  df[is.na(df)] <- 0
  df
}

ESD <- function(fit, baseline = "pooled", std = FALSE) {
  MNIeffect(fit, baseline, std, measure = "ESD")
}

RESD <- function(fit, baseline = "pooled", std = FALSE) {
  MNIeffect(fit, baseline, std, measure = "RESD")
}