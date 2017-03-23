# Filename: WaveNumL.R
# 
# Author: Luke Miller  Mar 22, 2017
###############################################################################

#' WaveNumL
#' 
#' A function to calculate wavenumber.
#' 
#' @param f A numeric vector of wave frequencies
#' @param h A numeric vector of water depths (usually in units of meters)
#' @return The wave number. 
#' @references Original MATLAB function by Urs Neumeier 
#' http://neumeier.perso.ch/matlab/waves.html
#' @author George Voulgaris, SUDO, 1992

WaveNumL <- function(f, h) {
	w <- 2*pi*f
	dum1 <- (w^2)*h/9.81
	dum2 <- dum1 + 
			(1.0+0.6522*dum1 + 0.4622*dum1^2 + 0.0864*dum1^4 + 
				0.0675*dum1^5)^(-1);
	dum3 <- sqrt(9.81*h*dum2^(-1)) / f
	y <- 2*pi*dum3^(-1)
	y # return y
}
