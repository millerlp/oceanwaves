---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


`oceanwaves` provides a set of functions to calculate summary
statistics
for ocean waves, using a record of surface heights as input. For sea surface
heights
deri ved from bottom-mounted pressure transducers, the package also
contains a
function `prCorr` to correct for depth attenuation of the pressure signal. 

`waveStatsSP` produces wave height and period statistics using spectral analysis
methods, while `waveStatsZC` calculates additional wave height and period
statistics based on a zero-crossing algorithm. 

Pressure corrections and wave statistics functions were adapted from Urs
Neumeier's `waves` functions for MATLAB.
[http://neumeier.perso.ch/matlab/waves.html](http://neumeier.perso.ch/matlab/waves.html)  

To install this package from within R, first install the package `devtools`
[https://CRAN.R-project.org/package=devtools]{https://CRAN.R-project.org/package=devtools}
and then install this package from Github:


```r
install.packages('devtools')
library(devtools)
install_github('millerlp/oceanwaves')
```

