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
#' @param method A character string indicating which spectral analysis method
#' should be used. Choose one of "welchPSD" (default) or "spec.pgram".
#' @param kernel An object of class 'tskernel' that defines a smoother for use
#' with spec.pgram method. 
#' @param segments Numeric value indicating the number of windowing segments to
#' use with welchPSD method.
#' @param ... Additional arguments to be passed to spectral analysis functions, 
#' such as windowfun
#' @return Data frame of wave parameters based on spectral methods

waveStatsSP <- function(PT, Fs, method = c('welchPSD','spec.pgram'), 
		kernel = NULL, segments = NULL, ...){
	# Inputs
	# PT = vector of pressure transducer values, converted to seawater height
	#		(height of sea surface above pressure transducer, units = meters)
	# 
	method = match.arg(method, choices = c('welchPSD','spec.pgram'))
	
	min_frequency = 0.05;	#	% mininum frequency, below which no correction is applied (0.05)
	max_frequency = 0.33;		#% maximum frequency, above which no correction is applied (0.33)	
	
	#
	# Prepare data for spectral analysis
	
	m <- length(PT);		# Get length of surface height record

	#####################
	h <- mean(PT);			#% mean water depth
	
	# Call oceanwaves::detrendHeight() 
	# Function that computes the least-squares fit of a straight line to 
	# the data and subtracts the resulting function from the data. The resulting
	# values represent deviations of surface height from the mean surface 
	# height in the time series.
	PT <- oceanwaves::detrendHeight(PT); 

	# Convert timeseries of surface heights to a timeseries object
	xt = ts(PT, frequency = Fs)
	
	if (method == 'spec.pgram'){
		if (is.null(kernel)){
			kernelval = kernel('daniell',c(9,9,9)) # Set default
		} else if (class(kernel) == 'tskernel') {
			kernelval = kernel # Use the user's kernel values
		} else {
			stop("Kernel for spec.pgram must be of class 'tskernel'")
		}
		# Use spec.pgram to estimate power spectral density, with smoothers
		pgram <- spec.pgram(xt, kernel = kernelval, taper=0.1, plot=FALSE)
		pgramm0 = (Fs*mean(pgram$spec))
		deltaf = pgram$freq[2]-pgram$freq[1] # bandwidth (Hz)
		integmin=min(which(pgram$freq >= 0)); # this influences Hm0 and other wave parameters
		integmax=max(which(pgram$freq <= max_frequency*1.5 ));
		moment = vector(length = 7)
		# Calculate moments of the spectrum, from -2nd to 0th to 4th
		# For a spectrum, the 0th moment represents the variance of the data, and
		# should be close to the value produced by simply using var(PT)
		for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
			# Note that the pgram$spec values are multiplied by 2 to normalize them
			# in the same fashion as a raw power spectral density estimator
			moment[i+3]=sum(pgram$freq[integmin:integmax]^i*
							(2*pgram$spec[integmin:integmax]))*deltaf;
		}	
		# Peak period, calculated from Frequency at maximum of spectrum 
		Tp = 1/pgram$freq[which.max(pgram$spec)] # units seconds
	} else if (method == 'welchPSD'){
		if (is.null(segments)){
			# Set default segment length for windowing
			Noseg = 4
		} else {
			Noseg = segments
		} 
		M = 2 * (length(mywaves)/ (Noseg+1)) / Fs
		seglength =  M
		
		wpsd = welchPSD(xt, seglength = M,two.sided = FALSE,method='mean',
				windowingPsdCorrection=TRUE)
		# Remove the zero-frequency entry
		wpsd$frequency = wpsd$frequency[-1]
		# Remove the zero-frequency entry from power as well
		wpsd$power = wpsd$power[-1] 
		
		deltaf = wpsd$frequency[2]-wpsd$frequency[1] # delta-frequency
		integmin=min(which(wpsd$frequency >= 0)); # this influences Hm0 and other wave parameters
		integmax=max(which(wpsd$frequency <= max_frequency*1.5 ));
		moment = vector(length = 7)
		# Calculate moments of the spectrum, from -2nd to 0th to 4th
		# For a spectrum, the 0th moment represents the variance of the data, and
		# should be close to the value produced by simply using var(PT)
		for (i in seq(-2,4,by=1)) { # calculation of moments of spectrum
			# Note that the wpsd$power values are multiplied by 2 to normalize them
			# in the same fashion as a raw power spectral density estimator
			moment[i+3]=sum(wpsd$frequency[integmin:integmax]^i*
							(wpsd$power[integmin:integmax]))*deltaf;
		}
		# Peak period, calculated from Frequency at maximum of spectrum 
		Tp = 1/wpsd$frequency[which.max(wpsd$power)]  # units seconds
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

	# spectral width parameters
	EPS2 = (moment[0+3]*moment[2+3]/moment[1+3]^2-1)^0.5;  
	EPS4 = (1 - moment[2+3]^2/(moment[0+3]*moment[4+3]))^0.5;

    
	results = data.frame(h = h, Hm0 = Hm0, Tp = Tp, m0 = m0, T_0_1 = T_0_1,
			T_0_2 = T_0_2, EPS2 = EPS2, EPS4 = EPS4)
		
	results # Return data frame
}

################################################################################

