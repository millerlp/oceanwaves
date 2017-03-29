#' PrCorr
#' 
#' A function to correct for depth attenuation of a water surface height 
#' pressure signal. 
#' 
#' @param pt A vector of water surface heights (units of meters usually).
#' @parma Fs Sampling frequency (units of Hz). Normally 4Hz for an OWHL logger.
#' @param zpt Height of the pressure sensor above the seabed (units of meters).
#' @param M Length of time series segments that will be used in the detrending
#' and attenuation correction operations. 512 samples is the default, should be
#' an even number.
#' @param CorrLim [min max] frequency for attenuation correction (Hz, optional, 
#' default [0.05 0.33]).
#' @return A vector of the depth-corrected surface heights (units of meters 
#' usually). Any original trend in the input data (such as tide change) is 
#' present in the output data. The returned surface height fluctuations will 
#' typically be more extreme than the raw input surface heights. 
#' @references Original MATLAB function by Urs Neumeier 
#' http://neumeier.perso.ch/matlab/waves.html Modified from the Pcorr3.m 
#' function written by T. Mason, SOC, January 1997
#' 
#' Each segment of pt will be linearly detrended, corrected for attenuation,
#' and the linear trend will be added back to the returned data.
#' 
#' @examples 
#' data(wavedata)
#' corrected = PrCorr(wavedata$SurfaceHeightRaw.m, Fs = 4, zpt = 0.1)
#' # Plot the results
#' plot(x = wavedata$DateTime, y = corrected, type = 'l', 
#'  ylab='Surface Height, m', xlab = 'Time')
#' lines(x = wavedata$DateTime, y = wavedata$SurfaceHeightRaw.m, col = 'red')
#' legend('topleft',legend=c('Corrected','Raw'),col=c('black','red'),lwd = 2)

PrCorr <- function(pt, Fs, zpt, M = 512, CorrLim = c(0.05, 0.33) ){
	###### 
	# Inputs
	# pt	- A vector of surface height values (meters) derived from the 
	#			original pressure sensor time series (uncorrected for depth
	#			attenuation effects)
	# Fs	- Sampling frequency (Hz). Normally 4 Hz for OWHL logger
	# zpt	- Height of pressure sensor above seabed (meters)
	# M		- Length of time series segments that will be used in the 
	#			detrending and attenuation correction operations. 512 samples
	#			is the default. Should be an even number.
	# Corr_lim - [min max] frequency for attenuation correction (Hz, 
	#                optional, default [0.05 0.33])
	library(signal) # for hanning() function
	
	# normally the maximum attenuation correction should not be higher than 5
	max_attenuation_correction <- 5
	# % higher value to process boat-waves for Jean Ellis
	# % max_attenuation_correction <- 20; 
	
	# mininum frequency, below which no correction is applied (0.05) = 20 s period
	min_frequency <- CorrLim[1]
	# maximum frequency, above which no correction is applied (0.33) = 3 s period
	max_frequency <- CorrLim[2]
	
	H_with_NaN <- pt	# make a copy of the input depth values
	notNA <- which(!is.nan(pt)) # get the indices of good data (not NA)
	pt<-pt[notNA] # make a reduced copy with only the good rows of pt
	
	# do_detrend=isempty(h); # little h is the mean water depth. MATLAB version
	do_detrend <- TRUE # Need to insert code to check for presence of  
	# an argument 'h', h and skip the detrending if it
	# was provided. LPM note: h is not currently an argument
	
	m <- length(pt)	# m = number of samples (all the NAs are excluded)
	
	Noverlap <- M/2	# length of overlap of segments. M is an argument to the
	# function and has a default of 512
	N <- (ceiling(m/M)) * M # length of array, will be 
							# zero-padded to nearest multiple of M
	
	f <- c(NA, seq(1,(M/2),by = 1) * Fs/M)  
	# f will be a vector of frequencies, with a NA in the 1st position. 
	
	
	x <- seq(1, M, by = 1) # vector of indices from 1 to M
	
	H <- vector(mode="numeric", length = N) # Make a vector of zeros, length N
	
	overlap_window <- signal::hanning(M); # R version of MATLAB's hann() function,
								# returns a L-point symmetric Hann window that
								# is used as an array of coefficients to 
								# combine overlapping segments
	
	# Mirror the overlap values in the 2nd half of overlap_window
	overlap_window[((M/2)+1):length(overlap_window)] <- 1 - overlap_window[1:(M/2)]
	
	# Step through segments of pt vector, detrend the data if needed
	for (q in seq(1, N-Noverlap,by = Noverlap)){
		o <- min(c(q+M-1, m))
		ptseg <- pt[q:o]
		seg_len <- length(ptseg)
		
		if (do_detrend){
			# Fit a simple linear regression through the segment of sea surface
			# height data and extract the intercept + slope
			trend <- coef(lm(ptseg ~  x[1:seg_len]))
			
			# Calculate a depth h at the midpoint of the segment of 
			# data using the regression coefficients stored in 'trend' 
			# (i.e. intercept and slope) 
			h <- trend[1] + ( trend[2] * ( (seg_len+1)/2 ) )  
			
			# Remove the linear trend from the segment of sea surface heights
			ptseg <- ptseg - (trend[1] + (trend[2]* x[1:seg_len])) 
			
			# Calculate the wave number for each frequency in f, using the
			# mean sea surface height h for this segment of data
			K <- WaveNumL(f,h);  # WaveNumL function defined separately
			
			# Correction factor of spectrum for pressure
			Kpt <- cosh(K*zpt) / cosh(K*h) 
			# Recall that zpt is supplied as an argument to the function, 
			# giving the height of the pressure transducer above the seabed.
			
			# Set values for frequencies outside the desired range to 
			# 1 so that no attenuation is applied 
			# (i.e. attenuate 0.05-0.33Hz only)
			Kpt[f < min_frequency | f > max_frequency] <- 1 	
			
			Kpt[Kpt < (1/max_attenuation_correction)] <- 1 / 
					max_attenuation_correction	
			#   % correction factor never higher than max_attenuation_correction
			#   % linear decrease of correction above f>max_frequency
			
			fb_max <- max(which(f <= max_frequency))
			
			# matlab fix() function rounds toward zero, R trunc() function 
			# should substitute
			# The min() function will return the smaller of 
			# fb_max+fix(length(K)/10) or length(K)
			mymin <- min(c(fb_max+trunc(length(K)/10), length(K)))
			fKlin <- seq(fb_max,mymin, by = 1) 
			
			Kpt[fKlin] <- ((seq(length(fKlin),1,by= -1)) * 
						(Kpt[fb_max]-1)/length(fKlin)) + 1
		
			Kpt[1] <- 1 # Overwrites that NA that's been hanging out in row 1
			
			# 2nd half of series is symmetric
			Kpt[seq(M,((M/2)+2),by=-1)] <- Kpt[2:(M/2)] 	
		}  # End of if(do_detrend) statement
		
		if (seg_len < M){
			ptseg[(seg_len+1):M] <- 0 # pad ptseg with zeros if its length is 
									  # shorter than M
		}# end of if(seg_len < M)
		
		P <- fft(ptseg) # Calculate spectrum
		Pcor <- P / Kpt # Apply correction factor
		# The R fft(x, inverse=TRUE) function returns unnormalized series, so it
		# is necessary to divide by length(x)
		# The Re() function keeps only the real part of the result
		Hseg <- Re(fft(Pcor, inverse = TRUE) / length(Pcor)) 

		# Keep only the part that is the same length as the  original segment.
		Hseg <- Hseg[1:seg_len] 
		
		if (do_detrend){
			# Add linear regression values back on to the detrended data so that
			# they are once again trended.
			Hseg <- Hseg + (trend[1] + trend[2]*x[1:seg_len])
			# Hseg should now contain sea surface height values (meters)
		}
		
		# Insert the corrected sea surface heights back into the vector H
		# using the coefficients in overlap_window to taper the ends of the
		# the segment Hseg
		H[q:o] <- H[q:o] + Hseg * overlap_window[1:seg_len]
		
		# Handle the 1st segment in H 
		if (q == 1){
			H[1: min(c(Noverlap,seg_len))] <- Hseg[1: min(c(Noverlap,seg_len))]
		}
		
		# Deal with the tail end of H
		if ( ( (q+M) >= N) & (seg_len > Noverlap)){
			H[(q+Noverlap):o] <- Hseg[(Noverlap+1):length(Hseg)]
		}
		
	}   # End of for (q in seq(1,N-Noverlap,by=Noverlap)) loop
		
	H <- H[1:m] # truncate H back to original length m
	
	# Place the new corrected heights back into the
	# the positions in the original vector that had
	# good data. 
	H_with_NaN[notNA] <- H; 
	# And rename that to be H
	H <- H_with_NaN 
	
	H # return H with all of the depth-corrected surface heights and the 
	# missing values (NA) re-inserted. Should be the same length as the input
	# vector pt. Units should be the same as the original input values
	# in the vector pt (usually meters).
	
}  # end of PrCorr function

