# The R-Package `MNIeffect`

This package implemets the effect size measures of first-order Measurement Non-Invariance (MNI) in [Schuhbeck et al. (2026)](https://doi.org/10.31234/osf.io/sjtp3_v2) for linear common factor models with (jointly) normally distributed factors. It can be installed using the `remotes` package:

```R
# install.packages("remotes")
remotes::install_github("TiziSchuh/MNIeffect")
```

### Main Functions

* `EAD`: group and across-group MNI effect sizes on item/sumscore level
* `MNIdir`: corresponding MNI effect directions
* `ED`: the impact of MNI on observed item/sumscore means 

Effects are calculated w.r.t. a global baseline (see, [Schuhbeck et al., 2026](https://doi.org/10.31234/osf.io/sjtp3_v2)) by default.

## Code-Example: Holzinger-Swineford 

All functions in this package require a multigroup `lavaan` fit object to calculate MNI effect sizes.

```R
# library("lavaan")

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

# estimate multigroup model 
fit <- cfa(model, data = HolzingerSwineford1939, group = "school")
```

_Note: The following tables were modified to display values rounded to two decimal places for simplicity._

### Group & Across-Group Effects of MNI

You can report _Cohen's-d_-like effect sizes of MNI using `EAD` with standardization enabled. 

```R
# library("MNIeffect")
 
> EAD(fit, std = TRUE)
              X1   X2   X3   X4   X5   X6   X7   X8   X9 Sumscore
Pasteur     0.00 0.14 0.23 0.21 0.28 0.25 0.23 0.04 0.05     0.07
Grant-White 0.01 0.12 0.23 0.23 0.30 0.27 0.25 0.04 0.05     0.07
Total       0.00 0.13 0.23 0.22 0.29 0.26 0.24 0.04 0.05     0.07
```
 
`0.00` indicates measurement invariance; higher values indicate higher non-invariance.

### MNI Effect Directions

You can report corresponding effect directions using `MNIdir`.

```R
# library("MNIeffect")
 
> MNIdir(fit)
               X1    X2    X3    X4    X5    X6    X7    X8    X9 Sumscore
Pasteur      1.00 -0.64  0.94 -1.00 -1.00 -1.00  1.00  0.93  0.90    -0.88
Grant-White -1.00  0.78 -0.98  1.00  1.00  1.00 -1.00 -0.89 -0.85     0.94
Total        0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00     0.00
```

`+1.00` indicates systematically higher responses than the baseline.  
` 0.00` indicates that the non-invariance cancels out completely w.r.t. the baseline.  
`-1.00` indicates systematically lower responses than the baseline.

### The Impact of MNI on Observed Means

The impact of MNI on observed means with respect to the baseline, in this case the group Paster, can be calculated using `ED`.  

```R
# library("MNIeffect")
 
> ED(fit, baseline = "Pasteur")
               X1   X2    X3   X4   X5   X6    X7    X8    X9 Sumscore
Pasteur      0.00 0.00  0.00 0.00 0.00 0.00  0.00  0.00  0.00     0.00
Grant-White -0.01 0.22 -0.49 0.49 0.72 0.55 -0.51 -0.07 -0.09     0.79
Total       -0.01 0.10 -0.24 0.24 0.35 0.26 -0.25 -0.04 -0.04     0.38
```

For instance: The Grant-White sumscore value `0.79` indicates that the (sumscore-)scale mean of Grant-White students is 0.79 points higher compared to a scenario in which they were to follow the baseline response behavior, i.e., in which the scale were to measure the latent trait exactly as it does for Pasteur students.
In this sense, +0.79 points may be attributed to MNI in a scale-mean comparison x̄<sub>Grant-White</sub>-x̄<sub>Pasteur</sub>.

# Source Code for Schuhbeck et al. (2026)

In addition, his repository contains the source codefor the empirical illustration  (see `illustration/`)  in:

Schuhbeck, T. M. B., Partsch, M. V., Sterner, P., & Goretzko, D. (2026, August 6). Effect-Based Assessments of Intergroup Comparability – The Practical Impact of Measurement Non-Invariance & Differential Item Functioning on Psychological Tests and Scales. Retrieved from osf.io/preprints/psyarxiv/sjtp3_v2
