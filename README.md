# `MNIeffect` Package

This package implements the effect size measures of first-order Measurement Non-Invariance (MNI) in [Schuhbeck et al. (2026)](https://doi.org/10.31234/osf.io/sjtp3_v2) for linear common factor models with (jointly) normally distributed factors. Installation using `remotes`:

```R
# install.packages("remotes")
remotes::install_github("TiziSchuh/MNIeffect")
```

_Note: Bug reports, suggestions and feedback are always welcome -> schuhbeck@psych.uni-frankfurt.de_

### Main Functions

* `EAD`: group and across-group MNI effect sizes on item/sumscore level
* `MNIdir`: corresponding MNI effect directions
* `ED`: the impact of MNI on observed item/sumscore means 

Effects are calculated w.r.t. a global baseline (see, [Schuhbeck et al., 2026](https://doi.org/10.31234/osf.io/sjtp3_v2)) by default.


## [^1] Code-Example: Holzinger-Swineford


[^1]: The tables in the code example were modified to display values rounded to two decimal places for simplicity.


All functions in this package require a multigroup `lavaan` fit object to calculate MNI effect sizes.

```R
# library("lavaan")
# library("MNIeffect")

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

# estimate multigroup model 
fit <- cfa(model, data = HolzingerSwineford1939, group = "school")
```

### Group & Across-Group Effects of MNI

You can report _Cohen's-d_-like effect sizes of MNI as `EAD` in the metric of global (i.e., across-group pooled) standard deviations.

```R
> EAD(fit, std = "global")

              X1   X2   X3   X4   X5   X6   X7   X8   X9 Sumscore
Pasteur     0.00 0.14 0.23 0.21 0.28 0.25 0.23 0.04 0.05     0.07
Grant-White 0.01 0.12 0.23 0.23 0.30 0.27 0.25 0.04 0.05     0.07
Total       0.00 0.13 0.23 0.22 0.29 0.26 0.24 0.04 0.05     0.07
```
 
`0.00` indicates measurement invariance; higher values indicate higher non-invariance.

### MNI Effect Directions

You can report corresponding effect directions using `MNIdir`.

```R 
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
> ED(fit, baseline = "Pasteur")

               X1   X2    X3   X4   X5   X6    X7    X8    X9 Sumscore
Pasteur      0.00 0.00  0.00 0.00 0.00 0.00  0.00  0.00  0.00     0.00
Grant-White -0.01 0.22 -0.49 0.49 0.72 0.55 -0.51 -0.07 -0.09     0.79
Total       -0.01 0.10 -0.24 0.24 0.35 0.26 -0.25 -0.04 -0.04     0.38
```

For instance: The Grant-White sumscore value `0.79` indicates that the (sumscore-)scale mean of Grant-White students is 0.79 points higher compared to a scenario in which they were to follow the baseline response behavior, i.e., in which the scale were to measure the latent trait exactly as it does for Pasteur students.
In this sense, +0.79 points may be attributed to MNI in a scale-mean comparison x̄<sub>Grant-White</sub>-x̄<sub>Pasteur</sub>.

# References

Schuhbeck, T. M. B., Partsch, M. V., Sterner, P., & Goretzko, D. (2026, November 14). Effect-Based Assessments of Intergroup Comparability – The Practical Impact of Measurement Non-Invariance & Differential Item Functioning on Psychological Tests and Scales. Retrieved from osf.io/preprints/psyarxiv/sjtp3_v2 [^2]

Schuhbeck, T. M. B., Sterner, P., & Goretzko, D. (2026). Quantifying Measurement Non-Invariance Beyond Simple Structure: The Closed Formulas of Universal Effect Size Measures for MI. Structural Equation Modeling: A Multidisciplinary Journal, 33(1), 81-88.

[^2]: The source code for the empirical illustration can be found in this repository under `illustration/`.

