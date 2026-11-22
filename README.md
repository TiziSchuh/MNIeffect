# The R-Package `MNIeffect`

This package implemets the effect size measures of Schuhbeck et al. (2026) – see the [Preprint on OSF](https://doi.org/10.31234/osf.io/sjtp3_v2) – 
for  linear factor models with (jointly) normally distributed factors. Installation of the package:

```R
# install.packages("remotes")
remotes::install_github("TiziSchuh/MNIeffect")
```


### Holzinger-Swineford Example 

A table of standardized group (EAD) and across-group (TEAD) effect sizes of measurement non-invariance on both item- and sumscore-level,

|            |   X1|   X2|   X3|   X4|   X5|   X6|   X7|   X8|   X9| Sumscore|
|:-----------|----:|----:|----:|----:|----:|----:|----:|----:|----:|--------:|
|Pasteur     | 0.00| 0.14| 0.23| 0.21| 0.28| 0.25| 0.23| 0.04| 0.05|     0.07|
|Grant-White | 0.01| 0.12| 0.23| 0.23| 0.30| 0.27| 0.25| 0.04| 0.05|     0.07|
|Total       | 0.00| 0.13| 0.23| 0.22| 0.29| 0.26| 0.24| 0.04| 0.05|     0.07|

obtained by the code

```R
# library("lavaan")
# library("MNIeffect")

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

# multigroup model 
fit <- cfa(model, data = HolzingerSwineford1939, group = "school")

# dataframe of rounded non-invariance effects 
df <- round(
        EAD(fit, std = "global"),
        digits = 2)
df
```

# Source Code for Schuhbeck et al. (2026)

In addition, his repository contains the source codefor the empirical illustration  (see `illustration/`)  in:

Schuhbeck, T. M. B., Partsch, M. V., Sterner, P., & Goretzko, D. (2026, August 6). Effect-Based Assessments of Intergroup Comparability – The Practical Impact of Measurement Non-Invariance & Differential Item Functioning on Psychological Tests and Scales. Retrieved from osf.io/preprints/psyarxiv/sjtp3_v2
