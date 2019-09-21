# Filename: joinOWHLfiles.R
# 
# Author: Luke Miller  Nov 28, 2018
###############################################################################


##########################
# Function joinOWHLfiles

#' Concatenate multiple OWHL raw csv data files
#' 
#'  Concatenate multiple OWHL raw csv data files  
#' 
#' Concatenate multiple OWHL raw csv data files into a single data frame and 
#' interpolate any missing 1-second intervals. The resulting data frame can 
#' be used to generate wave statistics from other functions in the package. 
#' 
#' The input csv files can have mission info in the first row, and then the 
#' following columns: POSIXt, DateTime, frac.seconds, Pressure.mbar, TempC
#' 
#' 
#' @param filenames A vector of OWHL csv filenames, including path.
#' @param timezone Specifies the time zone the raw data timestamps represent. UTC is 
#' the preferred time zone for simplicity.
#' @param verbose Logical argument TRUE or FALSE, specifying if progress messages should be output. 
#' 
#' @return A data frame with the individual csv files concatenate together
#' 
#' @export
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar


joinOWHLfiles <- function(filenames, timezone = 'UTC', verbose = TRUE) {
	# Get a list of all of the csv files in the directory
	# filelist <- dir(filedir, pattern = '*.csv', full.names=TRUE)
		
	# There is a little bit of mission info stored in the first row
	missioninfo <- scan(filenames[1], what = character(),
			nlines = 1, sep = ',')
	
	if (verbose) {
	  message("Loading data files")
	  pb <- txtProgressBar(min = 0, max = length(filenames), style = 3)
	}
	# Open the raw data files and concatenate them.
	for (f in 1:length(filenames)){
		if (verbose){ setTxtProgressBar(pb,f) }
		# To get the real column headers, skip the first line
		dattemp <- read.csv(filenames[f], skip = 1)
		###########################
# Columns:
# POSIXt: elapsed seconds since 1970-01-01 00:00:00 (unix epoch) in whatever
#         timezone the sensor was set to during deployment
# DateTime: human-readable character date and time, in whatever timezone the
#         sensor was set to during deployment (Pacific Daylight Time for the
#         example file 15032600_halfday.CSV)
# frac.seconds: integer value representing the fractional seconds of the 
#         sample * 100. To get the actual fractional second, divide these values
#         by 100 (resulting in 0.0, 0.25, 0.50, 0.75)
# Pressure.mbar: temperature-compensated pressure in millibar, from the 
#         MS5803-14BA sensor. This is absolute pressure (including atmospheric)
# TempC:  internal temperature of the MS5803-14BA pressure sensor.
		#########################
# Convert the DateTime column to a POSIXct object. 

		dattemp$DateTime <- as.POSIXct(dattemp$DateTime, tz = timezone) 	
# Add on the fractional seconds for each row using the values in frac.seconds
		dattemp$DateTime <- dattemp$DateTime + (dattemp$frac.seconds/100)
		
		# Concatenate the data files. 
		if (f == 1){
			dat <- dattemp
		} else if (f > 1){
			dat <- rbind(dat,dattemp)
		}
	}
	
	if (verbose) { close(pb) } # shut off progress bar
	
# Reorder the concatenated data frame by the DateTime values in case the files
# were not fed in in chronological order
	dat <- dat[order(dat$DateTime),]
	
	###########
# Filter suspect time points. Suspect time points have ms values that aren't
# equal to 0, 25, 50, or 75.
# Make a matrix to hold the test results
	filt <- matrix(0,nrow = nrow(dat), ncol = 5)
	filt[,1] <- rep(0L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 0
	filt[,2] <- rep(25L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 25
	filt[,3] <- rep(50L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 50
	filt[,4] <- rep(75L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 75
# Sum every row of filt. If none of the tests above returned true, the result
# in the 5th column will be 0 (== FALSE). Otherwise it will be TRUE
	filt[,5] <- rowSums(filt[,1:4]) 
# Now find every row that has a suspect ms value
	badrows <- which(!(filt[,5]))
	if (length(badrows)>0){
		# Remove any rows with suspect ms values
		dat <- dat[-badrows,]	
		message('Removed ',length(badrows),' suspect rows from data.')
	} else {
		# do nothing if there were no bad rows
	}
	################################################################################
# If there are brief breaks in the data (1 second) due to SD card write issues
# then it would be nice to interpolate those to fill them and be able to use
# the whole chunk of time. 
	bounds <- which(diff(dat$DateTime) > 0.25 ) # Find all breaks > 0.25 seconds
	bounds2 <- which(diff(dat$DateTime) > 1.25) # Find all breaks > 1.25 seconds
	
# Define a function to find entries in bounds that aren't in bounds2
	"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
# Get a vector of indices in dat where the next time step is missing 1 second
	shortbreaks <- bounds %w/o% bounds2
	
	# If there are missing time steps of 1 second, proceed with filling them
	if (length(shortbreaks) > 0) {
	  # Determine the length of the new data frame with the added interpolated rows
	  newlength <- nrow(dat) + (length(shortbreaks) * 4)
	  # Allocate the new data frame
	  newdat <- data.frame(
	    POSIXt = numeric(length = newlength),
	    DateTime = numeric(length = newlength),
	    frac.seconds = integer(length = newlength),
	    Pressure.mbar = numeric(length = newlength),
	    TempC = numeric(length = newlength)
	  )
	  newdat$DateTime <-
	    as.POSIXct(newdat$DateTime, origin = '1970-1-1',
	               tz = timezone)
	  
	  nextEmptyRow <- 0 # counter variable
	  
	  
	  message("Interpolating missing seconds")
	  if (verbose) {
	    pb <- txtProgressBar(min = 0,
	                        max = length(shortbreaks),
	                        style = 3)
	  }
	  
	  
	  # Go through each location and insert the missing second of data (4 rows) via
	  # linear interpolation of the pressure and temperature values
	  for (i in 1:length(shortbreaks)) {
	    if (verbose) {
	      setTxtProgressBar(pb, i)
	    }
	    if (i == 1) {
	      # Grab the first chunk of good rows
	      tempdat <- dat[1:shortbreaks[i], ]
	      # Get the next 4 good rows after the skip
	      afterskip <- dat[(shortbreaks[i] + 1):(shortbreaks[i] + 4), ]
	      # Get record of last good pressure + temp, and next good pressure + temp
	      beforePress <- tempdat$Pressure.mbar[nrow(tempdat)]
	      afterPress <- afterskip$Pressure.mbar[1]
	      beforeTempC <- tempdat$TempC[nrow(tempdat)]
	      afterTempC <- afterskip$TempC[1]
	      # Fit a linear model to the pressure values
	      presslm <- coef(lm(c(beforePress, afterPress) ~ c(1, 6)))
	      xs <- seq(2, 5)
	      # Interpolate the 4 missing values
	      interpPress <- (xs * presslm[2]) + presslm[1]
	      # Do the same routine for temperature
	      templm <- coef(lm(c(beforeTempC, afterTempC) ~ c(1, 6)))
	      interpTempC <- (xs * templm[2]) + templm[1]
	      
	      afterskip2 <- afterskip
	      afterskip2$POSIXt <- afterskip2$POSIXt - 1
	      afterskip2$DateTime <- afterskip2$DateTime - 1
	      afterskip2$Pressure.mbar <- interpPress
	      afterskip2$TempC <- interpTempC
	      
	      # Add the missing 4 rows to newdat
	      newdat[1:(shortbreaks[1]), ] <-
	        tempdat  # First copy over the good data
	      newdat[(shortbreaks[1] + 1):(shortbreaks[1] + 4),] <- afterskip2
	      # Update the nextEmptyRow to the first empty row in newdat
	      nextEmptyRow <- shortbreaks[1] + 5
	      
	      # Move on to the 2nd entry in shortbreaks, which should get added to
	      # the newdat data frame
	    } else {
	      # Grab the next chunk of good rows
	      tempdat <- dat[(shortbreaks[i - 1] + 1):shortbreaks[i], ]
	      # Get the next 4 good rows after the skip
	      afterskip <- dat[(shortbreaks[i] + 1):(shortbreaks[i] + 4), ]
	      # Get record of last good pressure + temp, and next good pressure + temp
	      beforePress <- tempdat$Pressure.mbar[nrow(tempdat)]
	      afterPress <- afterskip$Pressure.mbar[1]
	      beforeTempC <- tempdat$TempC[nrow(tempdat)]
	      afterTempC <- afterskip$TempC[1]
	      # Fit a linear model to the pressure values
	      presslm <- coef(lm(c(beforePress, afterPress) ~ c(1, 6)))
	      xs <- seq(2, 5)
	      # Interpolate the 4 missing values
	      interpPress <- (xs * presslm[2]) + presslm[1]
	      # Do the same routine for temperature
	      templm <- coef(lm(c(beforeTempC, afterTempC) ~ c(1, 6)))
	      interpTempC <- (xs * templm[2]) + templm[1]
	      
	      afterskip2 <- afterskip
	      afterskip2$POSIXt <- afterskip2$POSIXt - 1
	      afterskip2$DateTime <- afterskip2$DateTime - 1
	      afterskip2$Pressure.mbar <- interpPress
	      afterskip2$TempC <- interpTempC
	      
	      # Add the good data to newdat. The value of nextEmptyRow currently
	      # should hold the 1st empty row of newdat, and we want to add on
	      # all of the good rows in tempdat.
	      newdat[nextEmptyRow:(nextEmptyRow + nrow(tempdat) - 1),] <-
	        tempdat
	      # Update nextEmptyRow now that new data have been inserted
	      nextEmptyRow <- nextEmptyRow + (nrow(tempdat))
	      # Add the interpolated rows in as well
	      newdat[(nextEmptyRow):(nextEmptyRow + 3), ] <- afterskip2
	      # Update the row counter to the current empty row of newdat
	      nextEmptyRow <- nextEmptyRow + 4
	      
	    }
	  }
	  
	  if (verbose) {
	    close(pb) # shut off progress bar
	    message('Filled ',i,' one second gaps')
	  }
	  # Finish the operation by appending the remaining good rows from dat
	  tempdat <- dat[(shortbreaks[i] + 1):nrow(dat), ]
	  newdat[nextEmptyRow:(nextEmptyRow + nrow(tempdat) - 1), ] <- tempdat
	} else {
	  # If there were no 1 second gaps, just return a copy of dat
	  newdat <- dat
	}
	# Return newdat
	newdat
}


