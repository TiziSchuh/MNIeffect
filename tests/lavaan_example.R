library("lavaan")
source("R/lavaan_wrapper.R")

# Holzinger Swineford
test_model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
 
fit <- cfa(test_model, data = HolzingerSwineford1939, group = "school")


# --- global baseline --- 

# via pooled parameters
round(ED(fit), digits = 2)
round(EAD(fit), digits = 2)
round(RESD(fit), digits = 2)

# via single group model (fit_global)
fit_global <- cfa(test_model, data = HolzingerSwineford1939, meanstructure = TRUE)
round(ED(fit, baseline = fit_global))
round(EAD(fit, baseline = fit_global))
round(RESD(fit, baseline = fit_global))

# --- reference group baseline ---

round(ED(fit, baseline = "Pasteur"), digits = 2)
round(EAD(fit, baseline = "Pasteur"), digits = 2)
round(RESD(fit, baseline = "Pasteur"), digits = 2)
