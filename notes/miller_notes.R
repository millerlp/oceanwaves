# TODO: implement windowed Welch's method-type spectral analysis
# TODO: fancy up the detrendHeight function to also return h (mean depth) and 
# coefficients so that I can eliminate the redundant code in prCorr.R
###########################################################################
setwd("D:/R_public/oceanwaves") # package must be the working directory
setwd("~/R_public/oceanwaves") # package must be the working directory
library(devtools)  # load devtools
devtools::document() # Regenerate documents, help files, namespace etc
build() # Generate a tar.gz file of the entire package
install() # Re-install oceanwaves after an update
load_all() # Actually load the functions in oceanwaves 

data(wavedata) # Load the example wave data

# Basic estimation of significant wave height using surface height variance
myx = detrendHeight(wavedata$SurfaceHeight.m)
varx = var(myx)  # Variance of surface height, after detrending
RawVarHsig = 4*sqrt(varx) # 4 times square root of surface height variance, sig H?

# Using the zero-crossing function to estimate sig. wave height. 
zcstats = waveStatsZC(wavedata$SurfaceHeight.m, Fs = 4)
periodZC = zcstats$Tsig
zerocrossHsig = zcstats$Hsig

# Spectral analysis method
library(astsa)
myspec = mvspec(wavedata$SurfaceHeight.m, 
                spans=c(13,11,3), taper = 0.5, plot = TRUE)
# Extract the frequency associated with the peak spectral density
myfreq = myspec$freq[which.max(myspec$spec)]
# Convert to a period
myperiod = 1/myfreq
# Since the data were collected at 4Hz, the period will be 4x higher than
# the real value, so divide by 4 to get the actual peak period in seconds
specPeriod = myperiod / 4

#####################
# Try welchPSD function from bspec package
#install.packages('bspec')
library(bspec)
data(lh)
tsval = lh
tsval = as.ts(wavedata$SurfaceHeight.m, start=c(1,1), frequency = 4)
spec1 = empiricalSpectrum(tsval,two.sided=TRUE)
plot(spec1$frequency[-1],spec1$power[-1],log='y',type='l')
spec2 = spectrum(tsval,plot=FALSE)
lines(spec2$freq,spec2$spec,col='red')
spec3 = welchPSD(tsval,seglength=450,method='mean')
lines(spec3$frequency,spec3$power, col = 'blue')


###############################################################################
# From https://cran.r-project.org/web/packages/psd/vignettes/normalization.pdf
# Normalization of Power Spectral Density estimates
set.seed(1234)
N <- 1028
x <- rnorm(N, mean = 0, sd = 1)
# Load signal processing
library(signal, warn.conflicts=FALSE)
# Construct an FIR filter
f <- c(0, 0.2, 0.2, 0.3, 0.3, 0.5)*2
m <- c(0, 0, 1, 1, 0, 0)
fw <- signal::fir2(N, f, m)
# complex filter response
fh <- signal::freqz(fw, Fs=1)
f <- c(0, 0.12, 0.12, 0.22, 0.22, 0.5)*2
fwl <- signal::fir2(N, f, m)
fhl <- signal::freqz(fwl, Fs=1)
f <- c(0, 0.28, 0.28, 0.38, 0.38, 0.5)*2
fwu <- signal::fir2(N, f, m)
fhu <- signal::freqz(fwu, Fs=1)
# convolution
xf <- signal::filter(fw, x)
# PSD using stats::spectrum
Sx <- spectrum(x, pad=1, plot=FALSE, taper=0.2)
Sxf <- spectrum(xf, pad=1, plot=FALSE, taper=0.2)
plot(Sx$freq,Sx$spec,type='l',col='grey70', log='y',ylim=c(1e-12,1e1))
lines(Sxf$freq,Sxf$spec,col='black')

xv = var(x)  # x should be from rnorm(length = 1028, mean=0, sd=1) above
X = fft(x) # carry out Fast Fourier Transform
# 
Sa = Mod(X) # amplitude spectrum (if you square this, you get energy spectral density)
Sp = Arg(X) # phase spectrum
XC = Conj(X) # conjugate of spectrum
all.equal(Se <- Sa^2, Se_2 <- Mod(XC * X), Se_2R <- Mod(X * XC))

fsamp <- 1 # sampling freq, e.g. Hz
fNyq <- fsamp/2 # Nyquist freq
Nf <- N/2 # number of freqs
nyfreqs <- seq.int(from=0, to=fNyq, length.out=Nf)
S <- Se[2:(Nf+1)] * 2 / N # Finally, the PSD!

# 0) Setup optimization function for dof, using conjugate gradients\\
# min L1 |PSD - Chi^2(dof)|
Chifit <- function(PSD){optim(list(df=mean(PSD)), function(dof){
				sum(log(PSD)) - sum(log(dchisq(PSD, dof))) }, method="CG")}
# 1) run optimization
Schi <- Chifit(S)
# Get 'df', the degrees of freedom, should be very close to the mean mSn below
print(dof <- Schi$par[[1]])
# compare with the mean and median
c(mSn <- mean(S), median(S))

mSn <- dof
test_norm <- function(sval, nyq, xvar){svar <- sval * nyq; return(svar/xvar)}
print(xv_1 <- test_norm(mSn, fNyq, xv))
## [1] 1.003214 Should be close to 1 if the two variance estimates match
xv_2 <- sum(S)/Nf * fNyq / xv # an alternate test
all.equal(xv_1, xv_2)


# Change sampling frequency and re-scale to get properly normalized spectra
fsamp <- 4  #20
fNyq <- fsamp / 2
freqs <- fsamp * nyfreqs # rescaling original frequencies from 1Hz to 20Hz
Snew <- S / fsamp # S was our original PSD from above
# Test variance crudely
mSn <- mean(Snew)
test_norm(mSn, fNyq, xv)
## [1] 0.9990345

fsamp <- 4  #20
xt <- ts(x, frequency=fsamp)
pgram20 <- spec.pgram(xt, pad=1, taper=0, plot=FALSE)
pgram01 <- spec.pgram(ts(xt, frequency=1), pad=1, taper=0, plot=FALSE)
mSn/mean(pgram20$spec) # This will be off by a factor of 2
pgram20Norm = (2*mean(pgram20$spec)) # Normalize spectrum() output by multiplying by 2
mSn/pgram20Norm # This should be ~1.0


################################################################################
# Wave data PSD
# Calculate variance of wave data (detrended to zero mean)
# The code below pretty much replicates the output of the matlab wavesp.m 
# spectral analysis methods, when the matlab spectrum options are set to 
# not window the data at all (bw = 0.001111). Note that the T_pc peak period
# will be unrealistic with this method, and appears to require at least some
# windowing to bring it into the realistic realm.
N = length(wavedata$SurfaceHeight.m) # sample size
fsamp = 4 # sampling rate (Hz)
x = detrendHeight(wavedata$SurfaceHeight.m)
#x <- rnorm(N, mean = 0, sd = 1)  # alternate simple test with known variance
# Calculate variance of the timeseries using the var() function
xv = var(x)
# Calculate fast fourier transform
X = fft(x)
Sa = Mod(X) # calculate amplitude spectrum 
Se = Sa^2 # convert to energy spectral density

fNyq = fsamp/2 # Nyquist frequency (half of sampling rate)
Nf = N/2 # Half of sample size
# Generate a set of frequencies from 0 to Nyquist frequency
nyfreqs = seq.int(from=0, to=fNyq, length.out=Nf) # original 
nyfreqs = seq.int(from=0, to=fNyq, length.out=Nf+1) # Modified: matches wavesp.m
nyfreqs = nyfreqs[-1] # Modified: Remove first value if you do Nf+1 above
# Calculate the PSD from the energy spectral density
S <- Se[2:(Nf+1)] * 2 / N # original

Snew = S / fsamp # Rescale spectrum by sampling frequency (Hz)
myvar = function(S,Nf,fNyq){ # Calculate variance
	sum(S)/ Nf * fNyq
}
print(fftvar <- myvar(Snew,Nf,fNyq))
xv
fftvar/xv # should equal ~1 if the spectral method is matching the var() output

print(Hs <- 4 * sqrt(fftvar)) # Calculate sig. wave height as 4*sqrt(variance)
print(Hs2 <- 4 * sqrt(xv)) # Calculate sig. wave height as 4*sqrt(variance)
# Frequency at maximum of spectrum (converted to period)
Tp = 1/nyfreqs[which.max(Snew)] 

#% frequency range over which the spectrum is integrated for calculation of the moments
integmin=min(which(nyfreqs >= 0)); # this influences Hm0 and other wave parameters
min_frequency = 0.05 # from wavesp.m 
max_frequency = 0.33 # from wavesp.m 
integmax=max(which(nyfreqs <= max_frequency*1.5 ));

df = nyfreqs[2]-nyfreqs[1] # bandwidth (Hz) (should be same as values from diff(nyfreqs)
moment = vector(length = 7)
for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
	moment[i+3]=sum(nyfreqs[integmin:integmax]^i*Snew[integmin:integmax])*df;
}						
		
m0 = moment[3];
Hm0 = 4 * sqrt(m0)
T_0_1 = moment[3]/moment[1+3] # average period m0/m1
T_0_2 = (moment[0+3]/moment[2+3])^0.5 # average period (m0/m2)^0.5
T_pc = moment[-2+3]*moment[1+3]/(moment[0+3]^2) # calculated peak period
EPS2 = (moment[0+3]*moment[2+3]/moment[1+3]^2-1)^0.5;  # spectral width parameter
EPS4 = (1 - moment[2+3]^2/(moment[0+3]*moment[4+3]))^0.5; # spectral width parameter

# Plot the basic spectrum
plot(nyfreqs[-1],Snew[-1], type = 'n',col='blue',log='xy',
		ylab = expression(m^2/Hz), xlab = 'Frequency, Hz')
abline(v = (1/Tp), lty = 2, col = 'grey50') # peak period
ax2 = axTicks(side = 3, log=TRUE,nintLog = 8)
axis(side=3,at=ax2, labels = 1/ax2)
mtext(side = 3, text = 'Period, s', line = 2.75)
lines(nyfreqs[-1],Snew[-1], col = 'blue')
####


###
# Use spec.pgram, with smoothers
xt = ts(x, frequency = fsamp)
pgram20 <- spec.pgram(xt, spans = c(9,9,9), taper=0.1, plot=FALSE)
pgram20m0 = (fsamp*mean(pgram20$spec))
df = pgram20$freq[2]-pgram20$freq[1] # bandwidth (Hz) (should be same as values from diff(nyfreqs)
integmin=min(which(pgram20$freq >= 0)); # this influences Hm0 and other wave parameters
min_frequency = 0.05 # from wavesp.m 
max_frequency = 0.33 # from wavesp.m 
integmax=max(which(pgram20$freq <= max_frequency*1.5 ));
moment = vector(length = 7)
for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
	# Note that the pgram20$spec values are multiplied by 2 to normalize them
	# in the same fashion as a raw power spectral density estimator
	moment[i+3]=sum(pgram20$freq[integmin:integmax]^i*
					(2*pgram20$spec[integmin:integmax]))*df;
}	
m0 = moment[3];
Hm0 = 4 * sqrt(m0)
T_0_1 = moment[3]/moment[1+3] # average period m0/m1
T_0_2 = (moment[0+3]/moment[2+3])^0.5 # average period (m0/m2)^0.5
T_pc = moment[-2+3]*moment[1+3]/(moment[0+3]^2) # calculated peak period
EPS2 = (moment[0+3]*moment[2+3]/moment[1+3]^2-1)^0.5;  # spectral width parameter
EPS4 = (1 - moment[2+3]^2/(moment[0+3]*moment[4+3]))^0.5;
Tp = 1/pgram20$freq[which.max(pgram20$spec)] 
# Plot the smoothed periodogram, being sure to multiply x2 to match raw PSD
lines(pgram20$freq,2*pgram20$spec, col = rgb(0,1,0))

# End of matlab/spectrum comparison section
####################################################################



################################################################################
# Try the welchPSD windowing function to compute Power Spectral Density
library(bspec)
# This function should approximate the old MATLAB spectrum function, which also
# used Welch's method of windowing data. It is not an exact match for the 
# wavesp.m MATLAB output, but it's close. 
wpsd = welchPSD(xt, seglength=225,two.sided=FALSE,method='mean',
		windowingPsdCorrection=TRUE)
# Remove the zero-frequency entry
wpsd$frequency = wpsd$frequency[-1]
# Remove the zero-frequency entry from power as well
wpsd$power = wpsd$power[-1] #  * 2/fsamp # If using two.sided=TRUE, don't multiply by 2/fsamp

df = wpsd$frequency[2]-wpsd$frequency[1] # delta-frequency (should be same as values from diff(nyfreqs)
integmin=min(which(wpsd$frequency >= 0)); # this influences Hm0 and other wave parameters
min_frequency = 0.05 # from wavesp.m 
max_frequency = 0.33 # from wavesp.m 
integmax=max(which(wpsd$frequency <= max_frequency*1.5 ));
moment = vector(length = 7)
for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
	# Note that the wpsd$power values are multiplied by 2 to normalize them
	# in the same fashion as a raw power spectral density estimator
	moment[i+3]=sum(wpsd$frequency[integmin:integmax]^i*
					(2*wpsd$power[integmin:integmax]))*df;
}
m0 = moment[3];
Hm0 = 4 * sqrt(m0)
T_0_1 = moment[3]/moment[1+3] # average period m0/m1
T_0_2 = (moment[0+3]/moment[2+3])^0.5 # average period (m0/m2)^0.5
T_pc = moment[-2+3]*moment[1+3]/(moment[0+3]^2) # calculated peak period
EPS2 = (moment[0+3]*moment[2+3]/moment[1+3]^2-1)^0.5;  # spectral width parameter
EPS4 = (1 - moment[2+3]^2/(moment[0+3]*moment[4+3]))^0.5;
Tp = 1/wpsd$frequency[which.max(wpsd$power)] 
# Calculate variance as the sum of power times 2, divided by number of frequencies 
# times the Nyquist frequency (which is equal to the width of each frequency band)
welchvar = (2*sum(wpsd$power))/(length(wpsd$frequency))*(fsamp/2)
# The variance above will not be exactly the same as the var(xt) of the original
# full-length time series, since the welch's estimate is based on averaging
# the spectra from several shorter windows of the original full-length timeseries
# But it ought to be close
welchvar/var(xt)

# Plot the smoothed periodogram, being sure to multiply x2 to match raw PSD
lines(wpsd$frequency,wpsd$power, col = rgb(1,0,1))




################################################################################
library(ggplot2)

swHeight = wavedata$SurfaceHeight.m
seg_len = length(swHeight)
secs = seq(0,seg_len, by = 0.25)  # 4Hz data = 0.25s steps
x = seq(1,seg_len,by = 1)
# Remove any linear trend (i.e. tide) from sea water height
swHeight = detrendHeight(swHeight)

# Calculate the dominant wave period. 
# Start with the spectrum periodogram function. The spans are just a 
# guess, and should be adjusted some day (they affect smoothing)
#	myspec = spec.pgram(swHeight2, spans=c(13,11), plot = FALSE)
myspec = mvspec(swHeight, spans=c(13,11,3), taper = 0.05, plot = FALSE)
spec.df <- data.frame(freq =myspec$freq, spec = myspec$spec)
names(spec.df)[2] = "taper 0.05"
spec.df[,'taper 0.1'] =  mvspec(swHeight, spans=c(13,11,3), taper = 0.1, plot = FALSE)$spec
spec.df[,'taper 0.2'] =  mvspec(swHeight, spans=c(13,11,3), taper = 0.2, plot = FALSE)$spec
spec.df[,'taper 0.5'] =  mvspec(swHeight, spans=c(13,11,3), taper = 0.5, plot = FALSE)$spec

spec.df <- reshape2:::melt.data.frame(data = spec.df, 
                                      id.vars = "freq", value.name = "spec", variable.name = "Taper")


# Extract the frequency associated with the peak spectral density
myfreq = myspec$freq[which.max(myspec$spec)]
# Convert to a period
myperiod = 1/myfreq
# Since the data were collected at 4Hz, the period will be 4x higher than
# the real value, so divide by 4 to get the actual peak period in seconds
myperiod = myperiod / 4
myperiod


minPeriod = 2 # units of seconds
maxPeriod = 20 # units of seconds
# Data are at 4Hz interval, so highest possible frequency is half that,
# meaning 2 Hz interval, 0.5 of 1 second.
g.period = rev(c(2,4,6,8,10,12,14,16,18))
g.labels = rev(c('2s','4','6','8','10','12','14','16','18')) 

# Convert the period values into frequencies. Multiply by 1/4 to convert from
# seconds to 4Hz
g.freqs = 1/g.period * (1/4)


ggplot(data = subset(spec.df, freq > 1/(maxPeriod*4) & freq < 1/(minPeriod*4) )) + 
  geom_path(aes(x = freq, y = spec, color = Taper)) + 
  scale_x_log10("Period", breaks = g.freqs, labels = g.labels) +
  #		scale_y_continuous(limits = ylims) +
  ylab('Spectrum') +
  labs(color='Taper') + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0))  
  # geom_vline(xintercept = buoyfreq)
# cat('Direction: ',buoyDir,'deg,\tperiod: ',buoyperiod,'s\n')
cat('Time: ',strftime(wavedata[1,'DateTime'],format='%Y-%m-%d %H:%M',tz='UTC'),'\n')


################################################################################
# NDBC wave spectra info, bandwidths, from http://www.ndbc.noaa.gov/wave.shtml 

NDBCspectrum = data.frame(Band = seq(1,47,by=1),
                          CenterFreqBand.Hz = 
                            c(0.020, seq(0.0325,0.0925, by = 0.005),
                              seq(0.1000,0.3500, by = 0.01), 
                              seq(0.3650,0.4850, by = 0.02)), 
                          Bandwidth.Hz = c(0.010, rep(0.005,13),
                                           rep(0.010, 26),
                                           rep(0.020, 7)))

NDBCspectrumOld = data.frame(Band = seq(1,33,by=1),
                             CenterFreqBand.Hz = seq(0.03,0.35, by = 0.01),
                             Bandwidth.Hz = rep(0.01,33))

################################################################################
################################################################################
# Using one of the above spectrum data frames, step through the available
# frequencies in your spectrum output, and average spectrum estimates for 
# the relevant frequencies in each band define in NDBCspectrum. 

swHeight = wavedata$SurfaceHeight.m
swHeight = detrendHeight(swHeight)
myspec = mvspec(swHeight, spans=c(13,11,3), taper = 0.5, plot = FALSE)
spec.df <- data.frame(freq =myspec$freq, spec = myspec$spec)

# Make a copy to hold output data
mybandavg = NDBCspectrum
mybandavg$spec = NA
for (i in 1:nrow(NDBCspectrum)){
  lowerFreq = NDBCspectrum$CenterFreqBand.Hz[i] -
    (NDBCspectrum$Bandwidth.Hz[i]/2)
  upperFreq = NDBCspectrum$CenterFreqBand.Hz[i] +
    (NDBCspectrum$Bandwidth.Hz[i]/2)
  # Find the indices in spec.df that match these frequencies
  indx1 = which.min(abs(lowerFreq - spec.df$freq))
  indx2 = which.min(abs(upperFreq - spec.df$freq))
  # Calculate the average of the spectra values for the range of frequencies
  mybandavg$spec[i] = mean(spec.df$spec[indx1:indx2]) 
}


# Calculate the variance (m_0) as the summation, over all relevant frequency
# bands, of the spectrum for each frequency multiplied by its bandwidth.
m0 = sum(mybandavg$spec * mybandavg$Bandwidth.Hz)
# Calculate significant wave height as 4 * sqrt(m0)
NDBCbandAvgHsig = 4 * sqrt(m0)
# Extract the frequency associated with the peak spectral density
myfreq = mybandavg$CenterFreqBand.Hz[which.max(mybandavg$spec)]
# Convert to a period
myperiod = 1/myfreq
# Since the data were collected at 4Hz, the period will be 4x higher than
# the real value, so divide by 4 to get the actual peak period in seconds
NDBCbandAvgTsig = myperiod / 4
#myperiod

# Calculate the variance for the 'raw' spectrum (not averaged)
lowerFreq = NDBCspectrum$CenterFreqBand.Hz[1] -
  (NDBCspectrum$Bandwidth.Hz[1]/2)
upperFreq = NDBCspectrum$CenterFreqBand.Hz[nrow(NDBCspectrum)] +
  (NDBCspectrum$Bandwidth.Hz[nrow(NDBCspectrum)]/2)
# Find the indices in spec.df that match these frequencies
indx1 = which.min(abs(lowerFreq - spec.df$freq))
indx2 = which.min(abs(upperFreq - spec.df$freq))
rawvar = sum(spec.df$spec[indx1:indx2]* (diff(spec.df$freq)[1]))
rawHsig = 4 * sqrt(rawvar)
# Extract the frequency associated with the peak spectral density
rawfreq = spec.df$freq[which.max(spec.df$spec)]
# Convert to a period
rawperiod = 1/rawfreq
# Since the data were collected at 4Hz, the period will be 4x higher than
# the real value, so divide by 4 to get the actual peak period in seconds
rawBandAvgperiod = rawperiod / 4

# Grab the closest matching buoy time point for this time interval
# First find the midpoint of the time series for this chunk
# midTime = mean(tempdat$DateTime)
# # Find the closest matching time value in buoydata
# buoyindx = which.min(abs(buoydata$DateTime - midTime))
# # Grab the CDIP Hs and peak period values
# cdipHsig = buoydata$Hs.m[buoyindx]
# cdipTp = buoydata$Tp.sec[buoyindx]


# Plot the smoothed spectrum and the averaged spectrum
titlepaste = paste('CDIP Hs: ',cdipHsig,'m\t avg Hs: ', format(myHsig,dig=3),'m\t',
                   'raw Hs: ',format(rawHsig,dig=3),'m\n',
                   'CDIP Tp: ', cdipTp,'s\t my Tp: ', format(myperiod,dig=3),'s\t raw Tp: ',
                   format(rawperiod,dig=3),'s\n',
                   'Time offset: ', format(midTime - 
                                             buoydata$DateTime[which.min(abs(buoydata$DateTime - midTime))],
                                           dig = 3),', ',midTime)
plot(mybandavg$CenterFreqBand.Hz,mybandavg$spec,type='b', 
     main = titlepaste)
lines(spec.df$freq,spec.df$spec,col='red',type = 'l')


##########################################################

cat('zero-cross T: ',periodZC,"\nraw spectrum T: ",specPeriod,
    "\nband-avg T: ",NDBCbandAvgTsig,'\nraw band-avg T: ',rawBandAvgperiod,'\n')

cat('raw variance Hsig: ',RawVarHsig,"\nzero-cross Hsig: ",zerocrossHsig, 
    "\nraw spec Hsig: ", rawHsig, "\nband-avg Hsig: ", NDBCbandAvgHsig,'\n')




################################################################################
### TEST ONLY - COMPARE TO MATLAB wavesp.m output
# A spectral analysis with no kernel smoothing and no tapering
swHeight = wavedata$SurfaceHeight.m
swHeight = detrendHeight(swHeight)
Fs = 4 # Sampling frequency
# swHeight = ts(swHeight, frequency = Fs)

N = length(swHeight)
# Sum of squared amplitude
sum(abs(swHeight)^2)  # 107.6604
# Mean squared amplitude
sum(swHeight^2)/N  # 0.02990566
# Fit an untapered, unwindowed spectrum
myspec = mvspec(swHeight, taper = 0, plot = FALSE, detrend = FALSE)
# Integrate spectrum by just adding all spec values together
sum(myspec$spec) # Roughly 1/2 of sum of squared amplitude above  53.83023
sum(myspec$spec)*2  # 107.6605

mean(swHeight) # should be roughly zero
mean(swHeight) + myspec$spec[floor(N/2)] + 2*sum(myspec$spec[2:floor(N/2)-1]) # 107.6604

myspecB = spectrum(swHeight, plot = FALSE, taper = 0, method='pgram')
library(psd)
myspecBnorm = normalize(myspecB, Fsamp = 4, src = 'double.sided')
sum(myspecBnorm$spec)  # 107.6605


s.fft = fft(swHeight)
1/N * sum(Mod(s.fft[1:floor(N/2)+1])^2)
a.PSD = 1/N * 2 * Mod(s.fft[1:(floor(N/2)+1)])^2
a.PSD[1] = a.PSD[1]/2
a.PSD[floor(N/2)+1] = a.PSD[floor(N/2)+1]/2
sum(a.PSD)  # 107.6604

spec.df <- data.frame(freq =myspec$freq, spec = myspec$spec)
# Attempt to Normalize spectrum
# spec.df$spec2[1:nrow(spec.df)] = spec.df$spec[1:nrow(spec.df)] * 2 / Fs
spec.df$spec2[1:nrow(spec.df)] = spec.df$spec[1:nrow(spec.df)] * 2 / Fs

min_frequency = 0.05;	#	% mininum frequency, below which no correction is applied (0.05)
max_frequency = 0.33;	
p = spec.df$spec2[2:nrow(spec.df)]
Snf = p
F = spec.df$freq[2:nrow(spec.df)]
#	% frequency range over which the spectrum is integrated for calculation of 
# the moments
integmin = which.min(F[F >= 0])	 #% this influences Hm0 and other wave parameters
integmax = which.max(F[F <= (max_frequency*1.5)])

df = F[1];	#	% bandwidth (Hz)
moment = numeric(length = 7)

# When i starts at -2, the first operation will be to take the square root
# of every frequency F (F^-2), multiply that by the spectrum at each frequency,
# sum all of those values, and multiply by the bandwidth (difference between
# frequency bands, df). This should be equivalent to m_0, the zeroeth moment
# of the spectrum (effectively the variance). Each subsequent value of i 
# produces the next moment (m_1, m_2, m_3 etc.). 

for (i in -2:4){	#		% calculation of moments of spectrum
  moment[i+3]=sum( (F[integmin:integmax]^i) * Snf[integmin:integmax])*df
  cat('moment ',i,' ',moment[i+3],'\n')
}
m0 = moment[0+3]

# Compare to matlab output from wavesp.m, using wavedata.csv file
# [res,names,spect] = wavesp(pt,[],4,[0.05 0.33],4,'a');
# m-2  =  3.73739
# m-1  =  0.280685
# m0  =  0.0281617
# m1  =  0.00345435
# m2  =  0.000527829
# m3  =  9.85453e-05
# m4  =  2.13276e-05

# Compare to matlab output from wavesp.m, using wavedata.csv file
# [res,names,spect] = wavesp(pt,[],4,[0.05 0.33],4,'a');
# 
# Columns 1 through 6 # Values are from spectral method
# 'h'        'Hm0'    'Tp'    'm0'    'T_0_1'    'T_0_2'
# 9.55132   0.67125  13.23529 0.02816 8.15253 7.30436
# Columns 7 through 11       # Values after EPS4 are from zero-cross routine
# 'T_pc'    'EPS2'    'EPS4'    'H_significant'    'H_mean'
# 16.27859   0.49570   0.73221    0.63949         0.41774   
# Columns 12 through 15  # All values from zero cross routine
# 'H_10'    'H_max'    'T_mean'    'T_s'
# 0.75663   0.94534   7.71982  10.78846



# rawvar = sum(spec.df$spec * (diff(spec.df$freq)[1]))
rawvar2 = sum(spec.df$freq^0 * spec.df$spec2) * diff(spec.df$freq)[1]
rawvar = sum(spec.df$spec)* diff(spec.df$freq)[1]
rawHsigNoFilt = 4 * sqrt(rawvar)
# Extract the frequency associated with the peak spectral density
rawfreq = spec.df$freq[which.max(spec.df$spec)]
# Convert to a period
rawperiod = 1/rawfreq
# Since the data were collected at 4Hz, the period will be 4x higher than
# the real value, so divide by 4 to get the actual peak period in seconds
rawperiodNoFilt = rawperiod / 4
# Data are at 4Hz interval, so highest possible frequency is half that,
# meaning 2 Hz interval, 0.5 of 1 second.
g.period = rev(c(2,4,6,8,10,12,14,16,18))
g.labels = rev(c('2s','4','6','8','10','12','14','16','18')) 
# Convert the period values into frequencies. Multiply by 1/4 to convert from
# seconds to 4Hz
g.freqs = 1/g.period * (1/4)
ggplot(data = subset(spec.df, freq > 1/(maxPeriod*4) & freq < 1/(minPeriod*4) )) + 
  geom_path(aes(x = freq, y = spec)) + 
  scale_x_log10("Period", breaks = g.freqs, labels = g.labels) +
  #		scale_y_continuous(limits = ylims) +
  ylab('Spectrum') +
  labs(color='Taper') + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0)) 

################################################################################


#have_matlab() # test that MATLAB was found

matlabwaves <- function(waves, Fs = 4) {
	library(matlabr)
	options(matlab.path = "C:\\Program Files (x86)\\MATLAB\\R2013a Student\bin")
# Convert surface height timeseries into character vector for MATLAB
	matwaves = rvec_to_matlab(waves)
# Generate a script file to call inside MATLAB
	code = c("cd('D:/MATLAB/work/waves');",
			paste0("PT = ",matwaves), 
			paste0("[res,names]=wavesp(PT,[],",Fs,",'az');"), 
			"cd('D:/R_public/oceanwaves/notes');",
			"save('test.txt', 'res', '-ascii','-tabs')")
	res = run_matlab_code(code, verbose=FALSE, wait = TRUE) # Run the code in MATLAB
	# Check that output file has been produced, pause otherwise
	watchdog = 1 
	while(!file.exists('D:/R_public/oceanwaves/notes/test.txt')){
		Sys.sleep(1)
		watchdog = watchdog + 1
		if (watchdog > 5) {cat('test.txt not found\n'); break}
	}
	if (file.exists('D:/R_public/oceanwaves/notes/test.txt')){
# Read in the output file produced by wavesp function in MATLAB
	output = read.table(file = 'D:/R_public/oceanwaves/notes/test.txt',
			sep = '\t',
			col.names=c( 'h','Hm0','Tp','m0','T_0_1','T_0_2','T_pc',
					'EPS2','EPS4','H_significant','H_mean','H_10',
					'H_max','T_mean','T_s'),
			colClasses = rep('numeric',15))
	file.remove('D:/R_public/oceanwaves/notes/test.txt')
	output
	} else {
		stop("Could not find test.txt output from MATLAB")
	}
}
################################################################################
mywaves = wavedata$SurfaceHeight.m

# Function to run full comparison of MATLAB and R methods, and plot results
runcomparison <- function(mywaves, Fs = 4, plot = FALSE){
# Run the matlab analysis
	matwaves = matlabwaves(mywaves, Fs = Fs)
# Run the R zero-cross version
	Rzc = waveStatsZC(mywaves,Fs = Fs)
# Run the R spectral analysis version
	Rsp = waveStatsSP(mywaves,Fs = Fs)
# Reformat results for easy comparison
	comparo = data.frame(Stat = colnames(matwaves), 
			MATLAB = as.vector(unlist(matwaves[1,])),
			R = c(as.vector(unlist(Rsp[1,])),as.vector(unlist(Rzc[1:6]))))
	# Calculate difference between MATLAB results and R results
	comparo$diff = comparo$MATLAB - comparo$R
	# Also calculate significant wave height directly using the variance
	# of the wave time series
	comparo[16,1] = NA
	comparo[16,2] = NA # This value isn't calculated by MATLAB function
	comparo[16,3] = 4 * sqrt(var(mywaves))
	comparo[16,4] = NA
	# Round off all values to a reasonable precision. Any very small differences
	# between the two methods will show up as 0. 
	comparo[,2:4] = round(comparo[,2:4],4)
	
	########## 
	if (plot) {
		# Plot the comparison results
		# Begin with wave height estimates
		par(mfrow = c(2,1),mar = c(5,5,1,1))
		cexval = 1.5
		cexaxis = 1.3
		mybg = 'grey70'
		plot(x = 0, y = 0,
				type = 'n',
				ylim = c(0,max(comparo[c(2,10:13),2:3])),
#			ylim = c(0,2.5),
				xlim = c(0.5,5.5),
				xaxt = 'n',
				las = 1,
				ylab = 'Wave height, m', xlab = '',
				cex.lab = 1.3, cex.axis = cexaxis)
		rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4],col=mybg)
		grid(nx = NA, ny=NULL, col = 'white', lty = 2)
		# Plot zero-cross Hmean
		points(x = c(0.9,1.1), y = c(comparo[11,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot spectral Hsig estimate
		points(x = c(1.9,2.1),y = c(comparo[2,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot zero-cross Hsig estimate
		points(x = c(2.9,3.1),y = c(comparo[10,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot zero-cross H10 = average height of highest 10% of waves
		points(x = c(3.9,4.1),y = c(comparo[12,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot zero-cross Hmax = highest overall wave
		points(x = c(4.9,5.1),y = c(comparo[13,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot point representing the Hsig calculated using timeseries variance
		points(x = 3.5, y = 4 * sqrt(var(mywaves)), pch = 20, col = 'black',
				cex = cexval)
		text(x = 3.5, y = 4 * sqrt(var(mywaves)), pos = 4, 
				labels = expression(H[sig]==4%*%sqrt(m[0])))
		axis(side = 1, at = c(1,2,3,4,5), 
				label = c('Hmean\nzero cross',
						'Hs\nspectral','Hs\nzero cross',
						'H10%\nzero cross',
						'Hmax\nzero cross'),
				line = 0, padj = 0.5, tick = FALSE)
		mtext(side = 1, text = 'Wave height', line = 4, cex = 1.3)
		legend('bottomright',legend=c('Matlab estimate','R estimate'), 
				col = c('blue','red'), cex = 1.2, pch = c(19,18), 
				box.col = mybg, bg = mybg)
		box()
		###############
		# Plot wave period stats
		par(mar = c(6,5,1,1))
		plot(x = 0, y = 0,
				type = 'n',
				ylim = c(0,max(comparo[c(3,5,6,14,15),2:3])),
				#			ylim = c(0,2.5),
				xlim = c(0.5,5.5),
				xaxt = 'n',
				las = 1,
				ylab = 'Wave period, s', xlab = '',
				cex.lab = 1.3, cex.axis = cexaxis)
		rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4],col='grey70')
		grid(nx = NA, ny=NULL, col = 'white', lty = 2)
		# Plot spectral T_0_1 estimate
		points(x = c(0.9,1.1), y = c(comparo[5,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot spectral T_0_2 estimate
		points(x = c(1.9,2.1),y = c(comparo[6,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot zero-cross mean period T_mean
		points(x = c(2.9,3.1),y = c(comparo[14,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot zero-cross T_s = period of highest 1/3 of waves
		points(x = c(3.9,4.1),y = c(comparo[15,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		# Plot spectral peak period Tp
		points(x = c(4.9,5.1),y = c(comparo[3,2:3]), col = c('blue','red'),
				pch = c(19,18), cex = cexval)
		axis(side = 1, at = c(1,2,3,4,5), 
				label = c('APD\nspectral\nm_0/m_1',
						'APD\nspectral\nsqrt(m0/m2)',
						'T_mean\nzero cross','T_sig\nzero cross',
						'T_peak\nspectral'),
				line = 0, padj = 0.5, tick = FALSE)
		mtext(side = 1, line = 4, text = 'Wave period', cex = 1.3)
		legend('bottomright',legend=c('Matlab estimate','R estimate'), 
				col = c('blue','red'), cex = 1.2, pch = c(19,18),
				 box.col = mybg, bg = mybg)
		box()
	}
	
	comparo # Return data frame
}
#runcomparison(wavedata$SurfaceHeight.m, Fs = 4, plot = TRUE)
################################################################################
# load Rdata file containing a data frame 'dat' with processed surface height
# time series from the Elsmore deployment 201610. 
load('D:/Dropbox/OWHL_misc/Deployment_Elsmore_201610.rda')
t1 = 1
# Go through chunks of data and compare MATLAB vs R wave stats results
runcomparison(dat$swDepthcorr.m[t1:(t1+3600)], Fs = 4, plot = TRUE);t1 = t1+3600

