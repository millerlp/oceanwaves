# Filename: wavesp.R
# 
# Author: Luke Miller  Mar 22, 2017
###############################################################################


###############################################################################
###############################################################################
## Function: wavesp 
## **********CURRENTLY NONFUNCTIONAL + UNTESTED****************

## Spectral wave parameters from pressure transducer data (PT) after 
## correction of the pressure attenuation (see function pr_corr)
#
## Ported to R from MATLAB by Luke Miller 2016
## % written by Urs Neumeier, 2003
## % modified from Coast07.m written by T Mason, UPCE, 1999 (psd modified 7/99)


wavesp <- function(PT, Zpt, Fs){
	# Inputs
	# PT = vector of pressure transducer values, converted to seawater height
	#		(height of sea surface above pressure transducer, units = meters)
	# 
	
	min_frequency = 0.05;	#	% mininum frequency, below which no correction is applied (0.05)
	max_frequency = 0.33;		#% maximum frequency, above which no correction is applied (0.33)	
	
	#%*************************************************************
	#		% Prepare data for spectral analysis
	
	m <- length(PT);		# Get length of surface height record
	# matlab fix() function rounds toward zero, R trunc() function 
	# should substitute
#	M = fix(m/Noseg/2)*2;	#% length of segment
	M <- trunc(m/Noseg/2)*2;	#% length of segment
	h <- mean(PT);			#% mean water depth
	PT <- detrend(PT); # MATLAB detrend() function computes the least-squares 
	# fit of a straight line (or composite line for piecewise linear trends) to 
	# the data and subtracts the resulting function from the data. The 
	# detrend() function defined separately in this source file does the
	# equivalent operation. 
	
	#%***********************************************************%
	#	% Calculation of spetral density 
	myspec = spec.pgram(PT, detrend = FALSE)
	p = myspec$spec
	Snf = p
	# Produce the vector of frequencies in Hz from the spectrum object
	F = myspec$freq * 4
	
	# The spectrum function in MATLAB is doing a routine where it steps through
	# the input data (PT) in windows of length M (900 in Noseg=4) and 
	# calculating the spectrum using overlapping windows (overlap =450, total 
	# of 7 windows in a 3600 sample chunk when M = 900). The question is whether
	# it is worthwhile to implement something similar in R
	
#	[P,F] <- spectrum(PT,M,M/2,[],Fs);  # % gives spectrum from 0 to Fnyq
	# The above in MATLAB feeds an empty vector [] to the window argument,
	# which results in the spectrum() function using hanning(M) to generate
	# the hanning window. But the hanning() function in MATLAB does not 
	# included the zero-weighted samples at the start and end of the output
	# window, while the hann() function does. 
	
#	p <- P(2:end,1)*2/Fs;	#% normalization so that psd=cov(7/99), and remove first						
#	F(1) <- [];		#% element of P and F, which contains no useful information
	
#	% frequency range over which the spectrum is integrated for calculation of 
# the moments
	integmin=which.min(F[F >= 0])	 #% this influences Hm0 and other wave parameters
	integmax=which.max(F[F <= (max_frequency*1.5 )])
	
	# TODO: I stopped here, since the matlab and R routines are deviating 
# quite a bit at this point, and numbers aren't coming up similar at all. Take
# a look back at the spectrum stuff above, since that's where they really 
# deviate.
	df = F[1];	#	% bandwidth (Hz)
	moment = numeric(length = 7)
	for (i in -2:4){	#		% calculation of moments of spectrum
		moment[i+3]=sum( (F[integmin:integmax]^i) * Snf[integmin:integmax])*df
	}
	#%	fprintf('m%g  =  %g\n',i,moment(i+3));
	
	m0 = moment[0+3]
	
}


