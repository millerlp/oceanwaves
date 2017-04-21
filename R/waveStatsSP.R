# Filename: waveStatsSP.R
# 
# Author: Luke Miller  2017-04-20, 2017
###############################################################################


###############################################################################
###############################################################################
## Function: waveStatsSP 

## Spectral wave parameters from pressure transducer data (PT) after 
## correction of the pressure attenuation (see function pr_corr)
#
## Ported to R from MATLAB by Luke Miller 2017
## % written by Urs Neumeier, 2003
## % modified from Coast07.m written by T Mason, UPCE, 1999 (psd modified 7/99)

#' Calculate ocean wave parameters using spectral analysis
#' 
#' Calculate ocean wave parameters using spectral analysis
#' 
#' Carries out spectral analysis of ocean wave height time series to estimate
#' common wave height statistics, including peak period, average period, 
#' and significant wave height.
#' 
#' @param PT A vector of surface heights that constitute a time series of 
#' observations. Units = meters
#' @param Fs Sampling frequency of the surface heights data. Units = Hz, i.e.
#' samples per second.
#' @return Data frame of wave parameters based on spectral methods

waveStatsSP <- function(PT, Fs){
	# Inputs
	# PT = vector of pressure transducer values, converted to seawater height
	#		(height of sea surface above pressure transducer, units = meters)
	# 
	
	min_frequency = 0.05;	#	% mininum frequency, below which no correction is applied (0.05)
	max_frequency = 0.33;		#% maximum frequency, above which no correction is applied (0.33)	
	
	#
	# Prepare data for spectral analysis
	
	m <- length(PT);		# Get length of surface height record
	####################
	# TODO: implement windowing method, such as Welch's method. 
	# matlab fix() function rounds toward zero, R trunc() function 
	# should substitute
#	M = fix(m/Noseg/2)*2;	#% length of segment
#	M <- trunc(m/Noseg/2)*2;	#% length of segment
	#####################
	h <- mean(PT);			#% mean water depth
	
	# Call oceanwaves::detrendHeight() 
	# Function that computes the least-squares fit of a straight line to 
	# the data and subtracts the resulting function from the data. The resulting
	# values represent deviations of surface height from the mean surface 
	# height in the time series.
	PT <- oceanwaves::detrendHeight(PT); 
	
	
	# Use spec.pgram to estimate power spectral density, with smoothers
	# Convert timeseries of surface heights to a timeseries object
	xt = ts(PT, frequency = Fs)
	pgram <- spec.pgram(xt, spans = c(9,9,9), taper=0.1, plot=FALSE)
	pgramm0 = (Fs*mean(pgram$spec))
	df = pgram$freq[2]-pgram$freq[1] # bandwidth (Hz) (should be same as values from diff(nyfreqs)
	integmin=min(which(pgram$freq >= 0)); # this influences Hm0 and other wave parameters
	min_frequency = 0.05 # from wavesp.m 
	max_frequency = 0.33 # from wavesp.m 
	integmax=max(which(pgram$freq <= max_frequency*1.5 ));
	moment = vector(length = 7)
	# Calculate moments of the spectrum, from -2nd to 0th to 4th
	# For a spectrum, the 0th moment represents the variance of the data, and
	# should be close to the value produced by simply using var(PT)
	for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
		# Note that the pgram$spec values are multiplied by 2 to normalize them
		# in the same fashion as a raw power spectral density estimator
		moment[i+3]=sum(pgram$freq[integmin:integmax]^i*
						(2*pgram$spec[integmin:integmax]))*df;
	}	
	# Estimate variance of time series (moment 0)
	m0 = moment[3]; 
	# Estimate significant wave height based on spectral moment 0, units meters
	# This value is approximately equal to the average of the highest one-third
	# of waves in the time series.  
	Hm0 = 4 * sqrt(m0) 
	# T_0_1, average period m0/m1, units seconds. Follows National Data Buoy
	# Center's method for average period (APD)
	T_0_1 = moment[3]/moment[1+3] 
	# T_0_2, average period (m0/m2)^0.5, units seconds. Follows Scripp's 
	# Institute of Oceanography's method for calculating average period (APD)
	# for their buoys. 
	T_0_2 = (moment[0+3]/moment[2+3])^0.5
	# Calculated peak period, units seconds. This method appears unreliable
	# for the R spectrum and should not be trusted
#	T_pc = moment[-2+3]*moment[1+3]/(moment[0+3]^2)
	T_pc = NA
	# spectral width parameters
	EPS2 = (moment[0+3]*moment[2+3]/moment[1+3]^2-1)^0.5;  
	EPS4 = (1 - moment[2+3]^2/(moment[0+3]*moment[4+3]))^0.5;
	# Peak period, calculated from Frequency at maximum of spectrum 
	Tp = 1/pgram$freq[which.max(pgram$spec)] # units seconds
    
	results = data.frame(h = h, Hm0 = Hm0, Tp = Tp, m0 = m0, T_0_1 = T_0_1,
			T_0_2 = T_0_2, T_pc = T_pc, EPS2 = EPS2, EPS4 = EPS4)
		
	results # Return data frame
}


