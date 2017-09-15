# Filename: plotWaveSpectrum.R
# 
# Author: Luke Miller  Apr 21, 2017
###############################################################################

#' Plot a basic spectrum 
#' 
#' @param freqspec A data frame containing a column of frequencies 'freq' and a
#' column of spectral power values 'spec'
#' @param Fs Frequency of sampled surface heights, units of Hz
#' @export

plotWaveSpectrum <- function(freqspec, Fs){
	# Extract the frequency associated with the peak spectral density
	maxfreq <- freqspec$freq[which.max(freqspec$spec)]
	# Convert to a period
	myperiod <- 1/maxfreq
	# Also grab max power
	maxpower <- max(freqspec$spec)
	minPeriod <- 2 # units of seconds
	maxPeriod <- 20 # units of seconds
# Data are at 4Hz interval, so highest possible frequency is half that,
# meaning 2 Hz interval, 0.5 of 1 second.
	g.period <- rev(c(2,4,6,8,10,12,14,16,18))
	g.labels <- rev(c('2s','4','6','8','10','12','14','16','18')) 
	
# Convert the period values into frequencies. 
	g.freqs <- 1 / g.period 
	
	myplot <- ggplot2::ggplot(data = subset(freqspec, freq > 1 / (maxPeriod) & 
									freq < 1 / (minPeriod) )) + 
			geom_path(aes(x = freq, y = spec)) + 
			scale_x_log10("Period, s", breaks = g.freqs, labels = g.labels,
					sec.axis = sec_axis(trans = ~., breaks = g.freqs,
							labels = format(g.freqs, dig=2, 
									drop0trailing = TRUE), 
					name = "Frequency, Hz")) +
			ylab(expression(Spectral~density*","~m^2/Hz)) +
			theme(plot.title = element_text(hjust = 0))  +
			annotate("text", x = maxfreq, y = maxpower, 
					label = paste('Tp = ',round(myperiod,1),'s'), hjust = 0,
					vjust = 0)
	print(myplot)
}
