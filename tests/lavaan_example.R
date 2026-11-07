library("lavaan")
source("R/lavaan_wrapper.R")

# Holzinger Swineford
test_model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(test_model, data = HolzingerSwineford1939, group = "school")
fit_global <- cfa(test_model, data = HolzingerSwineford1939, meanstructure = TRUE)

# effect of non-invariance w.r.t global fit
EAD(fit, baseline = fit_global)

# effect of non-invariance w.r.t reference group
EAD(fit, baseline = "Pasteur")