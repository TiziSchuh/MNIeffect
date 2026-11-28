<img align="right" src="tests/MNIeffect_v2.svg" width="100">

# `MNIeffect` Package



This package implements the effect size measures of first-order Measurement Non-Invariance (MNI) in [Schuhbeck et al. (2026a)](https://doi.org/10.31234/osf.io/sjtp3_v2) for linear common factor models with (jointly) normally distributed factors. Installation using `remotes`:

```R
# install.packages("remotes")
remotes::install_github("TiziSchuh/MNIeffect")
```

_Note: Bug reports, suggestions and feedback are always welcome -> schuhbeck@psych.uni-frankfurt.de_

# Quick Overview

All functions in this package require a multigroup `lavaan` fit object as their input and return a dataframe of non-invariance values relative to the chosen baseline (see documentation). 
The main functions are:
* `EAD` returns (_Cohen's-d_-like) group and (_Cohen's-f_-like) across-group MNI effect sizes on item/sumscore level
* `MNIdir` returns corresponding MNI effect directions
* `ED` returns impacts of MNI on observed item/sumscore means


## Code Example: Holzinger-Swineford

```R
# library("lavaan")
# library("MNIeffect")

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

# estimate multigroup model 
fit <- cfa(model, data = HolzingerSwineford1939, group = "school")
```

_Note: The following tables were modified to display values rounded to two decimal places for simplicity._

```R
> EAD(fit, std = "global")

              X1   X2   X3   X4   X5   X6   X7   X8   X9 Sumscore
Pasteur     0.00 0.14 0.23 0.21 0.28 0.25 0.23 0.04 0.05     0.07
Grant-White 0.01 0.12 0.23 0.23 0.30 0.27 0.25 0.04 0.05     0.07
Total       0.00 0.13 0.23 0.22 0.29 0.26 0.24 0.04 0.05     0.07
```
 
`0.00` indicates measurement invariance; higher values indicate higher non-invariance.


```R 
> MNIdir(fit)

               X1    X2    X3    X4    X5    X6    X7    X8    X9 Sumscore
Pasteur      1.00 -0.64  0.94 -1.00 -1.00 -1.00  1.00  0.93  0.90    -0.88
Grant-White -1.00  0.78 -0.98  1.00  1.00  1.00 -1.00 -0.89 -0.85     0.94
Total        0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00     0.00
```

`+1.00` indicates systematically higher responses than the baseline.  
` 0.00` indicates that the non-invariance (relative to the baseline) cancels out across latent levels.
`-1.00` indicates systematically lower responses than the baseline.


```R
> ED(fit, baseline = "Pasteur")

               X1   X2    X3   X4   X5   X6    X7    X8    X9 Sumscore
Pasteur      0.00 0.00  0.00 0.00 0.00 0.00  0.00  0.00  0.00     0.00
Grant-White -0.01 0.22 -0.49 0.49 0.72 0.55 -0.51 -0.07 -0.09     0.79
Total       -0.01 0.10 -0.24 0.24 0.35 0.26 -0.25 -0.04 -0.04     0.38
```

For instance: The Grant-White sumscore value `0.79` indicates that the (sumscore-)scale mean of Grant-White students is 0.79 points higher compared to a scenario in which they were to follow the baseline response behavior, i.e., in which the scale were to measure the latent trait exactly as it does for Pasteur students.
In this sense, +0.79 points may be attributed to MNI in a scale-mean comparison x̄<sub>Grant-White</sub>-x̄<sub>Pasteur</sub>.


# Documentation

The package contains five functions: `ED`, `EAD`, `MNIdir` and (for taxonomic completeness) `ESD` and `RESD`.  
All functions follow the same argument syntax with mandatory argument `fit`, a multigroup `lavaan` fit object.

## Functions

```R 
ED(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE)
```
returns the group effect **Expected Difference (ED)** and the across-group effect **Total Expected Difference (TED)** for item/sumscore responses.


```R 
EAD(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE)
```
returns the group effect **Expected Absolute Difference (EAD)** and the across-group effect **Total Expected Absolute Difference (TEAD)** for item/sumscore responses.

```R 
ESD(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE)
```
returns the group effect **Expected Squared Difference (ESD)** and the across-group effect **Total Expected Squared Difference (TESD)** for item/sumscore responses.

```R 
RESD(fit, baseline = "pooled", std = FALSE, include.sumscore = TRUE)
```
returns the group effect **Root Expected Squared Difference (RESD)** and the across-group effect **Root Total Expected Squared Difference (RTESD)** for item/sumscore responses.

```R 
MNIdir(fit, baseline = "pooled", include.sumscore = TRUE)
```
returns the **direction of** group and across-group **non-invariance** for item/sumscore responses.

## Optional Argument `baseline`

| Option                       |                                                                                |
| ---------------------------- | ------------------------------------------------------------------------------ |
| `baseline = "pooled"`        | Defines a global baseline by pooling intercepts and loadings across all groups |
| `baseline = my_fit`          | Defines the baseline via a single-group `lavaan` fit object                    |
| `baseline = "my_group_name"` | Defines a reference group baseline                                             |
 

## Optional Argument `std`

| Option                       |                                                                                       |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| `std = FALSE`                | Return values in the orignial item/sumscore metric                                    |
| `std = "global"`             | Return values in the metric of global (i.e., across-group-pooled) standard deviations |
| `std = "group"`              | Return group-effects in the metric of the respective group's standard deviation       |
| `std = "reference"`          | Return group-effects in the metric of the reference group's standard deviation        |
| `std = "pooled"`             | Return group-effects in the metric of (group-reference-)pooled standard deviations    |

For options `group`, `reference` and `pooled`, all across-group effects are reported in the metric of globally pooled standard deviations.  
Option `pooled` may be suitable, for instance, if pairwise comparisons to the reference group constitute the quantities of interest.

## Optional Argument `include.sumscore`

| Option                       |                                                                                        |
| ---------------------------- | -------------------------------------------------------------------------------------- |
| `include.sumscore = TRUE`    | Return values for item AND sumscore responses                                          |
| `include.sumscore = FALSE`   | Return values for item responses only                                                  |

# References

Schuhbeck, T. M. B., Partsch, M. V., Sterner, P., & Goretzko, D. (2026a, November 14). Effect-Based Assessments of Intergroup Comparability – The Practical Impact of Measurement Non-Invariance & Differential Item Functioning on Psychological Tests and Scales. Retrieved from osf.io/preprints/psyarxiv/sjtp3_v2 [^2]

Schuhbeck, T. M. B., Sterner, P., & Goretzko, D. (2026b). Quantifying Measurement Non-Invariance Beyond Simple Structure: The Closed Formulas of Universal Effect Size Measures for MI. Structural Equation Modeling: A Multidisciplinary Journal, 33(1), 81-88.

[^2]: The source code for the empirical illustration can be found in this repository under `illustration/`.

