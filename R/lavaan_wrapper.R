#' @noRd
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


#' Ths function implements a lavaan wrapper for the effect size measures
#' of measurement non-invariance in:
#' 
#' Schuhbeck, T. M. B., Partsch, M. V., Sterner, P., & Goretzko, D. 
#' (2026, November 14). Effect-Based Assessments of Intergroup Comparability 
#' – The Practical Impact of Measurement Non-Invariance & Differential 
#' Item Functioning on Psychological Tests and Scales. 
#' Retrieved from osf.io/preprints/psyarxiv/sjtp3_v2
#' 
#' It provides the functionality to calculate the effects 
#' (S)ED, (S)EAD, (S)ESD, (S)RESD, (S)TED, (S)TEAD, (S)TESD, (S)RTESD 
#' and their respective effect directions in lavaan for linear common 
#' factor models with (jointly) normally distributed factors, based on the
#' formulaic derivation in:
#' 
#' Schuhbeck, T. M. B., Sterner, P., & Goretzko, D. (2026). Quantifying 
#' Measurement Non-Invariance Beyond Simple Structure: The Closed Formulas 
#' of Universal Effect Size Measures for MI. Structural Equation Modeling: 
#' A Multidisciplinary Journal, 33(1), 81-88.
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Option `"pooled"`: uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param std Option `FALSE`: reports ED in orginial metric. 
#' Option `"global"`: reports ED in global/globally pooled SDs. 
#' Option `"group"`: reports ED in each group's respective SD.
#' Option `"reference"`: reports ED in reference group's SD.
#' Option `"pooled"`: reports ED in (reference-focal-)pooled SD.
#' NOTE: for conceptual reasons, total effects are ALWAYS standarized using 
#' global SDs for ANY option apart from `FALSE`.
#' @param include.sumscore Determines whether to include sumscore-level 
#' effects; can be `TRUE` or `FALSE`. 
#'  
#' @returns A dataframe of item-level (and sumscore-level) group and 
#' across-group non-invariance effects, either in the original metric 
#' or standardized
#' 
#' @noRd 
MNIeffect <- function(fit, baseline, std, include.sumscore, measure) {
  
  lambda <- lapply(lavInspect(fit, "est"), function(x) {unname(x$lambda)})
  nu <- lapply(lavInspect(fit, "est"), function(x) {unname(x$nu)})
  alpha <- lapply(lavInspect(fit, "est"), function(x) {unname(x$alpha)})
  psi <- lapply(lavInspect(fit, "est"), function(x) {unname(x$psi)})

  nitems <- length(nu[[1]])
  ngroups <- lavInspect(fit, "ngroups")
  nobs <- lavInspect(fit, "nobs")
  ntotal <- lavInspect(fit, "ntotal")
  p <- nobs/ntotal # group proportions

  if (inherits(baseline, "lavaan")) {
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
  

  get_sd <- function(sd.type) {
    # initialize with 1s = no standardization
    sd <- matrix(1, ngroups, nitems)
    sd.sumscore <- rep(1, ngroups)

    if (sd.type == "global") {
      if (inherits(baseline, "lavaan")) { 
        sd_base <- sqrt(diag(lavTech(baseline, "sampstat")[[1]]$cov))
        sd <- t(replicate(ngroups, sd_base))
        sd.sumscore <- replicate(ngroups,sd(get_sumscores(baseline)[[1]]))
      }
      else {
        sd_pool <- rep(0, nitems)
        sd.sumscore_pool <- 0 
        for (g in 1:ngroups) {
          sd_pool <- sd_pool + (nobs[g]-1) * diag(lavTech(fit, "sampstat")[[g]]$cov)
          sd.sumscore_pool <- sd.sumscore_pool + (nobs[g]-1) * (sd(get_sumscores(fit)[[g]])^2)
        }
        sd_pool <- sqrt(sd_pool/(sum(nobs)-ngroups))
        sd.sumscore_pool <- sqrt(sd.sumscore_pool/(sum(nobs)-ngroups))
        sd <- t(replicate(ngroups, sd_pool))
        sd.sumscore <- replicate(ngroups, sd.sumscore_pool)
      }
    }
    else if (sd.type == "group") {
      for (g in 1:ngroups) {
        sd[g,] <- sqrt(diag(lavTech(fit, "sampstat")[[g]]$cov))
        sd.sumscore[g] <- sd(get_sumscores(fit)[[g]])
      }
    }
    else if (sd.type == "reference") {
      sd_ref <- sqrt(diag(lavTech(fit, "sampstat")[[baseline]]$cov))
      sd <- t(replicate(ngroups, sd_ref))
      sd.sumscore_ref <- sd(get_sumscores(fit)[[baseline]])
      sd.sumscore <- replicate(ngroups, sd.sumscore_ref)
    }
    else if (sd.type == "pooled") {
      sd_ref <- sqrt(diag(lavTech(fit, "sampstat")[[baseline]]$cov))
      sd.sumscore_ref <- sd(get_sumscores(fit)[[baseline]])
      r <- match(baseline, lavInspect(fit, "group.label"))
      for (g in 1:ngroups) {
        sd[g,] <- (nobs[r]-1)*(sd_ref^2) + (nobs[g]-1)*diag(lavTech(fit, "sampstat")[[g]]$cov)
        sd[g,] <- sqrt(sd[g,]/(nobs[g]+nobs[r]-2))
        sd.sumscore[g] <- (nobs[r]-1)*(sd.sumscore_ref^2) + (nobs[g]-1)*(sd(get_sumscores(fit)[[g]])^2)
        sd.sumscore[g] <- sqrt(sd.sumscore[g]/(nobs[g]+nobs[r]-2))
      }
    }

    # for across-group effects
    sd_total <- matrix(1, nrow=1, nitems)
    sd.sumscore_total <- 1

    if (std != FALSE) {
      if (inherits(baseline, "lavaan")) { 
        sd_total <- sqrt(diag(lavTech(baseline, "sampstat")[[1]]$cov))
        sd.sumscore_total <- sd(get_sumscores(baseline)[[1]])
      }
      else {
        sd_pool <- rep(0, nitems)
        sd.sumscore_pool <- 0 
        for (g in 1:ngroups) {
          sd_pool <- sd_pool + (nobs[g]-1) * diag(lavTech(fit, "sampstat")[[g]]$cov)
          sd.sumscore_pool <- sd.sumscore_pool + (nobs[g]-1) * (sd(get_sumscores(fit)[[g]])^2)
        }
        sd_total <- sqrt(sd_pool/(sum(nobs)-ngroups))
        sd.sumscore_total <- sqrt(sd.sumscore_pool/(sum(nobs)-ngroups))
      }
    }

    list(
      group = sd,
      total = sd_total,
      group.sumscore = sd.sumscore,
      total.sumscore = sd.sumscore_total
    )
  } 

  sds <- get_sd(std) 
  sd <- sds$group
  sd_total <- sds$total
  sd.sumscore <- sds$group.sumscore
  sd.sumscore_total <- sds$total.sumscore


  get_effects <- function(effect.type) {
    # initialize with zeros
    effect <- matrix(0, ngroups, nitems)
    effect_total <- numeric(nitems)
    effect.sumscore <- matrix(0, ngroups, ncol = 1)
    effect.sumscore_total <- 0

    if(measure == "ED"){
      # item-level
      for (i in 1:nitems) {
        for (g in 1:ngroups) {
          effect[g,i] <- EDitem(nu_base[i], lambda_base[i, ], 
                                nu[[g]][i], lambda[[g]][i, ], 
                                alpha[[g]], psi[[g]]) 
          effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
          # standardize 
          effect[g,i] <- effect[g,i] / sd[g,i]
        }
        # standardize
        effect_total[i] <- effect_total[i] /sd_total[i]
      }
      # sumscore-level
      for (g in 1:ngroups) {
        effect.sumscore[g] <- EDtest(nu_base, lambda_base, 
                                    nu[[g]], lambda[[g]], 
                                    alpha[[g]], psi[[g]]) 
        effect.sumscore_total <- effect.sumscore_total + p[g] * effect.sumscore[g]
        # standardize
        effect.sumscore[g] <- effect.sumscore[g] / sd.sumscore[g] 
      }
      # standardize
      effect.sumscore_total <- effect.sumscore_total / sd.sumscore_total 
    }

    if(measure == "EAD"){ 
      # item-level
      for (i in 1:nitems) {
        for (g in 1:ngroups) {
          effect[g,i] <- EADitem(nu_base[i], lambda_base[i, ], 
                                nu[[g]][i], lambda[[g]][i, ], 
                                alpha[[g]], psi[[g]])
          effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
          # standardize 
          effect[g,i] <- effect[g,i] / sd[g,i]
        }
        # standardize
        effect_total[i] <- effect_total[i] /sd_total[i]
      }
      # sum-score level
      for (g in 1:ngroups) {
        effect.sumscore[g] <- EADtest(nu_base, lambda_base, 
                                    nu[[g]], lambda[[g]], 
                                    alpha[[g]], psi[[g]])
        effect.sumscore_total <- effect.sumscore_total + p[g] * effect.sumscore[g]
        # standardize
        effect.sumscore[g] <- effect.sumscore[g] / sd.sumscore[g] 
      }
      # standardize
      effect.sumscore_total <- effect.sumscore_total / sd.sumscore_total 
    }

    if(measure == "ESD"){ 
      # item-level
      for (i in 1:nitems) {
        for (g in 1:ngroups) {
          effect[g,i] <- ESDitem(nu_base[i], lambda_base[i, ], 
                                nu[[g]][i], lambda[[g]][i, ], 
                                alpha[[g]], psi[[g]])
          effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
          # standardize 
          effect[g,i] <- effect[g,i] / sd[g,i]
        }
        # standardize
        effect_total[i] <- effect_total[i] /sd_total[i]
      }
      # sumscore-level
      for (g in 1:ngroups) {
        effect.sumscore[g] <- ESDtest(nu_base, lambda_base, 
                                    nu[[g]], lambda[[g]], 
                                    alpha[[g]], psi[[g]])
        effect.sumscore_total <- effect.sumscore_total + p[g] * effect.sumscore[g]
        # standardize
        effect.sumscore[g] <- effect.sumscore[g] / sd.sumscore[g] 
      }
      # standardize
      effect.sumscore_total <- effect.sumscore_total / sd.sumscore_total 
    }

    if(measure == "RESD"){
      # item-level
      for (i in 1:nitems) {
        for (g in 1:ngroups) {
          effect[g,i] <- ESDitem(nu_base[i], lambda_base[i, ], 
                                nu[[g]][i], lambda[[g]][i, ], 
                                alpha[[g]], psi[[g]])
          effect_total[i] <- effect_total[i] + p[g] * effect[g,i]
          effect[g,i] <- sqrt(effect[g,i])
          # standardize 
          effect[g,i] <- effect[g,i] / sd[g,i]
        }
        effect_total[i] <- sqrt(effect_total[i])
        # standardize
        effect_total[i] <- effect_total[i] /sd_total[i]
      }
      # sumscore-level
      for (g in 1:ngroups) {
        effect.sumscore[g] <- ESDtest(nu_base, lambda_base, 
                                    nu[[g]], lambda[[g]], 
                                    alpha[[g]], psi[[g]])
        effect.sumscore_total <- effect.sumscore_total + p[g] * effect.sumscore[g]
        effect.sumscore[g] <- sqrt(effect.sumscore[g])
        # standardize
        effect.sumscore[g] <- effect.sumscore[g] / sd.sumscore[g] 
      }
      effect.sumscore_total <- sqrt(effect.sumscore_total)
      # standardize
      effect.sumscore_total <- effect.sumscore_total / sd.sumscore_total 
    }

    list(
      group = effect,
      total = effect_total,
      group.sumscore = effect.sumscore,
      total.sumscore = effect.sumscore_total
    )
  }

  effects <- get_effects(measure) 
  effect <- effects$group
  effect_total <- effects$total
  effect.sumscore <- effects$group.sumscore
  effect.sumscore_total <- effects$total.sumscore


  frame <- data.frame(
    rbind(effect, effect_total),
    row.names = append(lavInspect(fit, what = "group.label"), "Total")
  )
  if (include.sumscore == TRUE) {
    frame.sumscore <- data.frame(
      `Sumscore` = rbind(effect.sumscore, effect.sumscore_total)
    )
    frame <- cbind(frame,frame.sumscore)
  }

  frame
}



#' Average Difference of Expected Responses
#' 
#' This function calculates the Expected Difference (ED) and Total Expected
#' Difference (TED) of responses for linear common factor models with 
#' (jointly) normally distributed factors.  
#' The effect size DeltaMACS (Nye & Drasgow, 2011) is equivalent to
#' `ED(fit, baseline = "my_ref")`
#' in the `Sumscore` column 
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param std Default `FALSE`: report ED in orginial metric. 
#' Option `"global"`: report ED in global/globally pooled SDs. 
#' Option `"group"`: report ED in each group's respective SD.
#' Option `"reference"`: report ED in reference group's SD.
#' Option `"pooled"`: report ED in (reference-focal-)pooled SD.
#' NOTE: total effects are ALWAYS standarized using global SDs
#' for ANY option apart from `FALSE`!
#' @param include.sumscore Default `TRUE`. 
#' Additionally calculate group and across-group ED effects for sumscores.  
#' 
#' @examples
#' model <- ' visual  =~ x1 + x2 + x3
#'            textual =~ x4 + x5 + x6
#'            speed   =~ x7 + x8 + x9 '
#' dat <- lavaan::HolzingerSwineford1939
#' 
#' # multigroup model
#' fit <- lavaan::cfa(model, data = dat, group = "school")
#' 
#' # single-group model (only needed for baseline = fit_global)
#' fit_global <- lavaan::cfa(model, data = dat, meanstructure = TRUE)
#' 
#' ED(fit)
#' ED(fit, std = "global")
#' ED(fit, baseline = fit_global, std = "global") 
#' ED(fit, baseline = "Pasteur", std = "pooled", include.sumscore = FALSE)
#' 
#' @returns A dataframe of item-level and sumscore-level ED group effects  
#' and TED across-group effects in the original metric; for
#' standardized effects SED and STED, use the optional argument `std`.
#'  
#' @export
ED <- function(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE) {
  MNIeffect(fit, baseline, std, include.sumscore, measure = "ED")
}

#' Average Absolute Difference of Expected Responses
#' 
#' This function calculates the Expected Absolute Difference (EAD) and Total Expected
#' Absolute Difference (TEAD) of responses for linear common factor models with 
#' (jointly) normally distributed factors.  
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param std Default `FALSE`: report EAD in orginial metric. 
#' Option `"global"`: report EAD in global/globally pooled SDs. 
#' Option `"group"`: returns EAD in each group's respective SD.
#' Option `"reference"`: report EAD in reference group's SD.
#' Option `"pooled"`: report EAD in (reference-focal-)pooled SD.
#' NOTE: total effects are ALWAYS standarized using global SDs
#' for ANY option apart from `FALSE`!
#' @param include.sumscore Default `TRUE`. 
#' Additionally calculate group and across-group EAD effects for sumscores.  
#' 
#' @examples
#' model <- ' visual  =~ x1 + x2 + x3
#'            textual =~ x4 + x5 + x6
#'            speed   =~ x7 + x8 + x9 '
#' dat <- lavaan::HolzingerSwineford1939
#' 
#' # multigroup model
#' fit <- lavaan::cfa(model, data = dat, group = "school")
#' 
#' # single-group model (only needed for baseline = fit_global)
#' fit_global <- lavaan::cfa(model, data = dat, meanstructure = TRUE)
#' 
#' EAD(fit)
#' EAD(fit, std = "global")
#' EAD(fit, baseline = fit_global, std = "global") 
#' EAD(fit, baseline = "Pasteur", std = "pooled", include.sumscore = FALSE)
#' 
#' @returns A dataframe of item-level and sumscore-level EAD group effects  
#' and TEAD across-group effects in the original metric; for
#' standardized effects SEAD and STEAD, use the optional argument `std`.
#'  
#' @export
EAD <- function(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE) {
  MNIeffect(fit, baseline, std, include.sumscore, measure = "EAD")
}

#' The Direction of Measurement Non-Invariance
#' 
#' This function calculates the effect direction of measurement non-invariance 
#' for linear common factor models with (jointly) normally distributed factors.  
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param include.sumscore Default `TRUE`. 
#' Additionally calculate group and across-group directions for sumscore 
#' non-invariance.  
#' 
#' @examples
#' model <- ' visual  =~ x1 + x2 + x3
#'            textual =~ x4 + x5 + x6
#'            speed   =~ x7 + x8 + x9 '
#' dat <- lavaan::HolzingerSwineford1939
#' 
#' # multigroup model
#' fit <- lavaan::cfa(model, data = dat, group = "school")
#' 
#' # single-group model (only needed for baseline = fit_global)
#' fit_global <- lavaan::cfa(model, data = dat, meanstructure = TRUE)
#' 
#' MNIdir(fit)
#' MNIdir(fit, baseline = fit_global) 
#' MNIdir(fit, baseline = "Pasteur", include.sumscore = FALSE)
#' 
#' @returns A dataframe of item-level and sumscore-level group and 
#' across-group effect directions. 
#'  
#' @export
MNIdir <- function(fit, baseline = "pooled", include.sumscore = TRUE) {
  df <- ED(fit, baseline, std = FALSE, include.sumscore) / EAD(fit, baseline, std = FALSE, include.sumscore) 
  df[is.na(df)] <- 0
  df
}

#' Average Squared Difference of Expected Responses
#' 
#' This function calculates the Expected Squared Difference (ESD) and Total 
#' Expected Squared Difference (TESD) of responses for linear common factor 
#' models with (jointly) normally distributed factors. 
#' In almost all cases, RESD and RTESD should be calculated instead.
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param std Default `FALSE`: report ESD in orginial metric. 
#' Option `"global"`: report ESD in global/globally pooled SDs. 
#' Option `"group"`: returns ESD in each group's respective SD.
#' Option `"reference"`: report ESD in reference group's SD.
#' Option `"pooled"`: report ESD in (reference-focal-)pooled SD.
#' NOTE: total effects are ALWAYS standarized using global SDs
#' for ANY option apart from `FALSE`!
#' @param include.sumscore Default `TRUE`. 
#' Additionally calculate group and across-group ESD effects for sumscores.  
#' 
#' @examples
#' model <- ' visual  =~ x1 + x2 + x3
#'            textual =~ x4 + x5 + x6
#'            speed   =~ x7 + x8 + x9 '
#' dat <- lavaan::HolzingerSwineford1939
#' 
#' # multigroup model
#' fit <- lavaan::cfa(model, data = dat, group = "school")
#' 
#' # single-group model (only needed for baseline = fit_global)
#' fit_global <- lavaan::cfa(model, data = dat, meanstructure = TRUE)
#' 
#' ESD(fit)
#' ESD(fit, std = "global")
#' ESD(fit, baseline = fit_global, std = "global") 
#' ESD(fit, baseline = "Pasteur", std = "pooled", include.sumscore = FALSE)
#' 
#' @returns A dataframe of item-level and sumscore-level ESD group effects  
#' and TESD across-group effects in the original metric; for
#' standardized effects SESD and STESD, use the optional argument `std`.
#'  
#' @export
ESD <- function(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE) {
  MNIeffect(fit, baseline, std, include.sumscore, measure = "ESD")
}

#' Root Average Squared Difference of Expected Responses
#' 
#' This function calculates the Root Expected Squared Difference (RESD) and Root 
#' Total Expected Squared Difference (RTESD) of responses for linear common factor 
#' models with (jointly) normally distributed factors. 
#' The effect size dMACS (Nye & Drasgow, 2011) is equivalent to
#' `RESD(fit, baseline = "my_ref", std = "pooled", include.sumscore = FALSE)`
#' 
#' @param fit A lavaan multigroup fit object 
#' @param baseline Default `"pooled"` uses pooled intercepts and loadings 
#' as global baseline responses. A single group lavaan `fit` object can 
#' also be used to define global baseline responses. 
#' For baseline responses defined by a reference group, 
#' pass the name of the reference group as a string.
#' @param std Default `FALSE`: report RESD in orginial metric. 
#' Option `"global"`: report RESD in global/globally pooled SDs. 
#' Option `"group"`: returns RESD in each group's respective SD.
#' Option `"reference"`: report RESD in reference group's SD.
#' Option `"pooled"`: report RESD in (reference-focal-)pooled SD.
#' NOTE: total effects are ALWAYS standarized using global SDs
#' for ANY option apart from `FALSE`!
#' @param include.sumscore Default `TRUE`. 
#' Additionally calculate group and across-group RESD effects for sumscores.  
#' 
#' @examples
#' model <- ' visual  =~ x1 + x2 + x3
#'            textual =~ x4 + x5 + x6
#'            speed   =~ x7 + x8 + x9 '
#' dat <- lavaan::HolzingerSwineford1939
#' 
#' # multigroup model
#' fit <- lavaan::cfa(model, data = dat, group = "school")
#' 
#' # single-group model (only needed for baseline = fit_global)
#' fit_global <- lavaan::cfa(model, data = dat, meanstructure = TRUE)
#' 
#' RESD(fit)
#' RESD(fit, std = "global")
#' RESD(fit, baseline = fit_global, std = "global") 
#' RESD(fit, baseline = "Pasteur", std = "pooled", include.sumscore = FALSE)
#' 
#' @returns A dataframe of item-level and sumscore-level RESD group effects  
#' and RTESD across-group effects in the original metric; for
#' standardized effects SRESD and SRTESD, use the optional argument `std`.
#'  
#' @export
RESD <- function(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE) {
  MNIeffect(fit, baseline, std, include.sumscore, measure = "RESD")
}