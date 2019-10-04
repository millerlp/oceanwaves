## Resubmission

If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.

* I have added a citation in the DESCRIPTION Description field for a reference to a method used in the package. The reference is old (1979) and not archived online, so I cannot provide a doi link. Instead I have substituted the journal's ISSN number in the format suggested as "authors (year, ISBN:...)" in hopes that this is sufficient information.


Please make sure that you do not change the user's options, par or 
working directory. If you really have to do so, please ensure with an 
*immediate* call of on.exit() that the settings are reset when the 
function is exited, similar to this:
...
oldpar <- par(mfrow = c(1,2)) #line i
on.exit(par(oldpar)) #line i+1
...
f.i.: waveStatsZC.R

* I have added the suggested on.exit() call to waveStatsZC.R on the line immediately after the existing par() call. 



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
