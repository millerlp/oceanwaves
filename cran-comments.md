## Resubmission

found:
"function written by T. Mason"
"author George Voulgaris"
Please always add all authors and copyright holders in the Authors@R
field with the appropriate roles.

* I have modifed these attributions after contacting the listed
persons. These persons are now listed in the DESCRIPTION file
as contributors.

You write information messages to the console that cannot be easily
suppressed.
Instead of print()/cat() rather use message()/warning() or
if(verbose)cat(..) if you really have to write text to the console.
f.i.: joinOWHLfiles()

* I have modified joinOWHLfiles() to add a @param verbose to allow the
user to suppress txtProgressBar and messages. cat() functions have been converted to message(). 


Please always make sure to reset to user's options, wd or par after you
changed it in examples and vignettes.
F.i: vignette:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

* I have modified the function waveStatsZC() to reset user's par options after plotting. This should also correct the identified issue in the vignette.

## Test environments
* local OS X install, R 3.6.1
* local Windows 10 install, R 3.6.1 
* Windows via http://win-builder.r-project.org/ R Under development (unstable) (2019-09-20 r77199)
* Windows via http://win-builder.r-project.org/ R version 3.6.1 (2019-07-05)


## R CMD check results
NOTE: This is an initial submission. 
There were no other ERRORs, WARNINGs, or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.
