
setwd("~/R_public/oceanwaves") # package must be the working directory
library(devtools)  # load devtools
build()
install() # Re-install oceanwaves after an update
load_all() # Actually load the functions in oceanwaves 

data(wavedata) # Load the example wave data

# Basic estimation of significant wave height using surface height variance
myx = detrend(wavedata$SurfaceHeight.m)
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



################################################################################
library(ggplot2)

swHeight = wavedata$SurfaceHeight.m
seg_len = length(swHeight)
secs = seq(0,seg_len, by = 0.25)  # 4Hz data = 0.25s steps
x = seq(1,seg_len,by = 1)
# Remove any linear trend (i.e. tide) from sea water height
swHeight = detrend(swHeight)

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
swHeight = detrend(swHeight)
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
# A spectral analysis with no kernal smoothing and no tapering
swHeight = wavedata$SurfaceHeight.m
swHeight = detrend(swHeight)
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
