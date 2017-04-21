# Filename: waveStatsZC.R
# Based on the function zero_crossing.m by Urs Neumeier, v 1.06
# From http://neumeier.perso.ch/matlab/waves.html
# 
# 
# Author: Luke Miller  Mar 22, 2017
###############################################################################

#' Calculate wave statistics using zero-crossing method
#' 
#' Calculate ocean wave summary statistics, including significant
#' wave height and period. 
#' 
#' Based on an upward zero-crossing algorithm originally
#' provided by Urs Neumeier, v1.06. 
#' 
#' @param data A numeric vector of water surface height data. The data do not
#' need to be detrended prior to use.
#' @param Fs Sampling frequency of the data, in Hz. 
#' @param threshold The minimum height necessary for a zero-crossing event to 
#' be considered a wave. 
#' @param plot Set to TRUE if summary histograms of wave heights and wave
#' periods are desired.
#' @return A list object containing summary statistic values. Includes
#' significant wave height, mean wave height, 1/10th wave height, maximum wave
#' height, mean wave period, significant period.
#' @references Original MATLAB function by Urs Neumeier:  
#' http://neumeier.perso.ch/matlab/waves.html
#' @examples
#' data(wavedata)
#' waveStatsZC(data = wavedata$SurfaceHeight.m, Fs = 4)
#' @export

waveStatsZC <- function(data, Fs, threshold = NULL, plot = FALSE){
	# Sanity checks
	if (!is.vector(data)){
		error('data must be a vector of wave heights')
	}
	if (Fs < 0){
		error('Frequency Fs must be greater than 0')
	}
	
#% The function was written for zero upward-crossing. 
#% To have zero downward-crossing (recommended) leave the next line uncommented
	data=-data
	
	# detrendHeight function removes a linear trend and should result in a 
	# zero-mean timeseries
	data = oceanwaves::detrendHeight(data) 
	
	# Make a shorter version of data with 0's removed
	d0 = data[!almostZero(data,0)]	
	# Make a vector of indices the same length as 'data'
	back0 = 1:length(data)
	# Keep only indices that are non-zero in 'data'
	back0 = back0[!almostZero(data,0)]
 		
	# The code below multiplies each entry in d0 by the next entry in d0
	# and checks if the resulting value is negative (less than 0).
	# If two values are both positive heights, they will 
	# return a positive value, and if two values are both negative 
	# heights, then they will also return a positive value. But if a 
	# set of two values switches from positive to negative, or negative
	# to positive, the result will be a negative number (less than 0), and 
	# this must represent a crossing of the zero height level. 
	# The which() function then returns the indices where this is TRUE.
	f = which(d0[1:(length(d0)-1)] * d0[2:length(d0)] < 0)
	# Extract the indices in 'data' where the height transitioned across 0
	crossing = back0[f]
	
	# If the height record started off positive, then the first crossing must
	# be downward, and we will ignore it.
	if (data[1] > 0) {
		crossing = crossing[-1] # Remove the first crossing index
	}

	# Subset the crossing vector (containing indices) in steps of 2. The 
	# first entry should be a zero upward crossing, since we just removed the 
	# former first entry if it was a downward crossing. Every 2nd entry after 
	# that should also represent an upward zero crossing event. 
	crossing = crossing[seq(1,length(crossing), by = 2)]

	# Make a 4-column matrix to hold the height, crest, trough, and period data 
	# for each wave.
	wave = matrix(NA, nrow = length(crossing)-1, ncol = 4)
	
	for (i in 1:(length(crossing)-1)){
		# Get the highest surface height for the interval between each 
		# upward zero crossing
		wave[i,2] = max(data[crossing[i]:crossing[(i+1)]])
		# Get the lowest surface height for the interval between each 
		# upward zero crossing. This will be made positive for ease of use later
		wave[i,3] = -min(data[crossing[i]:crossing[(i+1)]])
	}

	# If no wave was found in data, then do nothing. Otherwise, process each
	# wave
	if (nrow(wave)>0){
		# Calculate the time between each upward zero crossing by calculating
		# the number of rows between each upward crossing and dividing
		# by the sampling rate (Fs, units of Hz). This will be 
		# the period for each wave.
		wave[,4] = diff(crossing)/Fs		
  
		if (is.null(threshold)){
			# If no minimum threshold for wave height was supplied, calculate a
			# threshold at 1% of the highest wave height (peak to trough) in the 
			# data set.
			threshold=0.01 * max(wave[,2]+wave[,3])
		} else if (threshold < 0){
			# If the supplied threshold was less than 0, throw an error
			stop('threshold must be greater than 0')
		}
		
		# Iterate through the 'wave' matrix and remove waves that are too small
		# and join them to the prior wave in the series
		 i = 0
		 while (i < nrow(wave)){
			 i = i + 1
			 # First check if the crest of wave i is less than threshold height
			# If the crest is too small, this will join the wave data with the
			# previous wave's data
			 if (wave[i,2] < threshold){
				 if (i != 1){
					# If we are on wave 2 or greater, rewrite wave i-1
					# with the max crest of wave i-1 or i, and the
					# max trough of wave i-i or i, and the
					# duration of wave i-1 and i.
					 wave[i-1,2:4] = c(max(wave[(i-1):i,2]),
							 max(wave[(i-1):i,3]),
							 sum(wave[(i-1):i,4]))
				 }
				 # Remove this row that was less than the threshold
				 wave = wave[-i,]
				 
			# Next check if the trough was smaller than the threshold
			# If the trough is small, then this will join the wave data with
			# the next wave in the series
			 } else if (wave[i,3] < threshold){
				 if (i != nrow(wave)){
					 # If we are on wave 2 or greater, rewrite wave i
					 # with the max crest of wave i or i+1, and the
					 # max trough of wave i or i+1, and the
					 # duration of wave i and i+1.
					 wave[i,2:4] = c(max(wave[i:(i+1),2]),
							 max(wave[i:(i+1),3]),
							 sum(wave[i:(i+1),4]))
					 # Remove the data for the next wave in the series
					 wave = wave[-(i+1),]
				 } else {
					 wave = wave[-i,]
				 }
			 }
		 }
		 # Add the crest and trough heights together for each wave to get
		 # overall wave height
		 wave[,1] = wave[,2] + wave[,3]
		
	} 	
					
# 'wave' has 1: wave height, 2: wave crest (Hcm), 3: wave trough (Htm), 4: period.

# Calculation of the wave statistics
	# Get the number of waves available to analyze
	nb=nrow(wave) 
	# Sort the wave matrix based on the 1st column (wave height), large to small
	wave = wave[order(wave[,1], decreasing = TRUE),]
	# Calculate the mean of the highest 1/3 of waves in the data set
	Hsig = mean(wave[1:round(nb*(1/3)), 1])
	# Overall mean wave height, for all waves (bigger than threshold)
	Hmean = mean(wave[ , 1])
	# Mean height of the upper 10% of all waves
	H10 = mean(wave[1:round(nb*0.1),1])
	# Maximum wave height
	Hmax = max(wave[,1])
	# Mean period
	Tmean = mean(wave[,4])
	# Mean period of Hsig (highest 1/3 of waves)
	Tsig = mean(wave[1:round(nb*(1/3)), 4])
	
	# Assemble a list to return
	waveStats = list(Hsig = Hsig, Hmean = Hmean, H10 = H10, Hmax = Hmax,
			Tmean = Tmean, Tsig = Tsig)
	
	if (plot){
		par(mfrow = c(2,1))
		hist(wave[,1], xlab = 'Wave height', main = '')
		hist(wave[,4], xlab = 'Wave Period', main = '')
		
	}

	# Return a list object with wave statistics
	waveStats

} # end of wave_stats_zc function



################################################################################ 
#' Test whether vector elements are effectively zero
#' 
#' Test whether elements in vector x are effectively zero,
#' within the square root of machine tolerance for a floating point number.
#' 
#' Returns TRUE for all vector elements that are within rounding error of zero.
#' 
#' @param x Numeric vector of values to compared against 0
#' @param y Value to be compared against. Default is 0. 
#' @param tolerance The maximum difference between x and y necessary to call a 
#' value non-zero
#' @return TRUE or FALSE for each element of the vector x. 

almostZero <- function(x, y = 0, tolerance=sqrt(.Machine$double.eps)) {
  diff <- abs(x - y)
  mag <- pmax(abs(x), abs(y))
  ifelse(mag > tolerance, diff/mag <= tolerance, diff <= tolerance)
}

