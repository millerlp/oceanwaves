
setwd("~/R_public/oceanwaves") # package must be the working directory
library(devtools)  # load devtools
build()
install() # Re-install oceanwaves after an update
load_all() # Actually load the functions in oceanwaves 

data(wavedata) # Load the example wave data

# Basic estimation of significant wave height using surface height variance
myx = detrend(wavedata$SurfaceHeight.m)
varx = var(myx)  # Variance of surface height, after detrending
4*sqrt(varx) # 4 times square root of surface height variance, sig H?

# Using the zero-crossing function to estimate sig. wave height. 
waveStatsZC(wavedata$SurfaceHeight.m, Fs = 4)
