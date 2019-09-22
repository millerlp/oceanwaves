
<!-- README.md is generated from README.Rmd. Please edit that file -->

`oceanwaves` provides a set of functions to calculate summary statistics
for ocean waves, using a record of surface heights as input. For sea
surface heights derived from bottom-mounted pressure transducers, the
package also contains a function `prCorr()` to correct for depth
attenuation of the pressure signal, and the `swDepth()` function from
the package `oce` can be used to convert pressure data into ocean
surface heights (see the included package vignette).

`waveStatsSP()` produces wave height and period statistics using
spectral analysis methods, while `waveStatsZC()` calculates additional
wave height and period statistics based on a zero-crossing algorithm.

An additional function `joinOWHLfiles` is provided to concatenate and
clean csv files from an Open Wave Height Logger pressure data logger, to
prepare the data for analysis with the functions described above. See
the package vignettes for example workflows to proceed from raw pressure
data csv files to summary wave statistics.

Pressure corrections and wave statistics functions were adapted from Urs
Neumeier’s `waves` functions for MATLAB, developed from earlier work by
Travis Mason and Magali Lecouturier.
<http://neumeier.perso.ch/matlab/waves.html>

To install this package from within R, first install the package
`devtools` <https://CRAN.R-project.org/package=devtools> and then
install this package from Github:

``` r
install.packages('devtools')
library(devtools)
install_github('millerlp/oceanwaves')
```
