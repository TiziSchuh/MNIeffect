library("lavaan")
source("R/lavaan_wrapper.R")

# Holzinger Swineford
test_model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

# multi-group model for group-specific response behavior 
fit <- cfa(test_model, data = HolzingerSwineford1939, group = "school")

# single-group model for global response behavior
fit_global <- cfa(test_model, data = HolzingerSwineford1939, meanstructure = TRUE)

# effect of non-invariance w.r.t global fit
ED(fit, baseline = fit_global)
EAD(fit, baseline = fit_global)
RESD(fit, baseline = fit_global)

# effect of non-invariance w.r.t reference group
EAD(fit, baseline = "Pasteur")