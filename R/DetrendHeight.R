# Filename: DetrendHeight.R
# 
# Author: Luke Miller  Mar 22, 2017
###############################################################################

#' Remove trend from a time series
#' 
#' Fits a straight line to a vector of values using lm(), and uses the
#' regression coefficients to subtract off the linear trend from the values.
#' Returns a series of residuals around the linear trend.
#' Typically this is used to remove a tidal trend from a ocean 
#' surface height time series before attempting to calculate statistics for 
#' waves, which are variations above and below the mean surface height. 
#' 
#' @param pt A vector of numeric values to be detrended
#' @return A vector of detrended values

DetrendHeight <- function(pt){
	# Determine length of pt
	seg_len <- length(pt)
	# Make a sequence of indices
	x <- seq(1,seg_len, by = 1)
	# Fit a simple linear regression through the segment of sea surface
	# height data and extract the intercept + slope
	trend <- coef(lm(pt ~  x[1:seg_len]))
	
	# Calculate a depth h at the midpoint of the segment of 
	# data using the regression coefficients stored in 'trend' 
	# (i.e. intercept and slope) 
	h <- trend[1] + ( trend[2] * ( (seg_len+1)/2 ) )  
	
	# Remove the linear trend from the segment of sea surface heights
	pt <- pt - (trend[1] + (trend[2]* x[1:seg_len])) 
	# Return the detrended vector
	pt
}
