library("lavaan")
source("R/lavaan_wrapper.R")

# Holzinger Swineford
test_model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
 
fit <- cfa(test_model, data = HolzingerSwineford1939, group = "school")
fit_global <- cfa(test_model, data = HolzingerSwineford1939, meanstructure = TRUE)

round(EAD(fit, std = "global"), digits = 3)
round(RESD(fit, std = "global"), digits = 3)
round(ESD(fit, baseline = "Pasteur", std = "global"), digits = 3)
round(RESD(fit, baseline = "Pasteur", std = "reference"), digits = 3)

MNIdir(fit)
round(MNIdir(fit, baseline = "Pasteur", include.sumscore = FALSE), digits = 3)

round(ED(fit, std = "global"), digits = 3)
round(ED(fit, std = "group"), digits = 3)
round(ED(fit, baseline = fit_global, std = "global"), digits = 3)
round(ED(fit, baseline = "Pasteur", std = "reference"), digits = 3)
round(ED(fit, baseline = "Pasteur", std = "pooled"), digits = 3)

round(EAD(fit), digits = 3)
round(EAD(fit, std = "group"), digits = 3)
round(EAD(fit, std = "global"), digits = 3)
round(EAD(fit, baseline = fit_global, std = "global"), digits = 3)
round(EAD(fit, baseline = "Pasteur", std = "pooled"), digits = 3)

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
