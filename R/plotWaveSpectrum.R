# Filename: plotWaveSpectrum.R
# 
# Author: Luke Miller  Apr 21, 2017
###############################################################################

#' Plot a basic smoothed spectrum 
#' 
#' @param PT A vector of surface heights
#' @param Fs Frequency of sampled surface heights, units of Hz

plotWaveSpectrum <- function(PT, Fs){
	# Use spec.pgram to estimate power spectral density, with smoothers
	# Convert timeseries of surface heights to a timeseries object
	xt = ts(PT, frequency = Fs)
	pgram <- spec.pgram(xt, spans = c(9,9,9), taper=0.1, plot=FALSE)
	
	pgram$spec = 2*pgram$spec # Multiply by 2 to match power spectral density
	myspec = data.frame(freq = pgram$freq,spec = pgram$spec)
	# Extract the frequency associated with the peak spectral density
	maxfreq = pgram$freq[which.max(pgram$spec)]
	# Convert to a period
	myperiod = 1/maxfreq
	# Also grab max power
	maxpower = max(pgram$spec)
	minPeriod = 2 # units of seconds
	maxPeriod = 20 # units of seconds
# Data are at 4Hz interval, so highest possible frequency is half that,
# meaning 2 Hz interval, 0.5 of 1 second.
	g.period = rev(c(2,4,6,8,10,12,14,16,18))
	g.labels = rev(c('2s','4','6','8','10','12','14','16','18')) 
	
# Convert the period values into frequencies. 
	g.freqs = 1/g.period 
	
	ggplot(data = subset(myspec, freq > 1/(maxPeriod) & freq < 1/(minPeriod) )) + 
			geom_path(aes(x = freq, y = spec)) + 
			scale_x_log10("Period", breaks = g.freqs, labels = g.labels) +
			#		scale_y_continuous(limits = ylims) +
			ylab(expression(Spectral~density*","~m^2/Hz)) +
			theme(plot.title=element_text(hjust=0))  +
			annotate("text", x = maxfreq, y = maxpower, 
					label = paste('Tp = ',round(myperiod,1),'s'), hjust = 0,
					vjust = 0)
}
