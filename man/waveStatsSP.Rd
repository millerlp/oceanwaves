% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waveStatsSP.R
\name{waveStatsSP}
\alias{waveStatsSP}
\title{Calculate ocean wave parameters using spectral analysis methods}
\usage{
waveStatsSP(
  data,
  Fs,
  method = c("welchPSD", "spec.pgram"),
  plot = FALSE,
  kernel = NULL,
  segments = NULL,
  ...
)
}
\arguments{
\item{data}{A vector of surface heights that constitute a time series of 
observations. Typical units = meters.}

\item{Fs}{Sampling frequency of the surface heights data. Units = Hz, i.e.
samples per second.}

\item{method}{A character string indicating which spectral analysis method
should be used. Choose one of \code{welchPSD} (default) or \code{spec.pgram}.}

\item{plot}{A logical value denoting whether to plot the spectrum. Defaults 
to \code{FALSE}.}

\item{kernel}{An object of class \code{tskernel} that defines a smoother for use
with \code{spec.pgram} method. If value is \code{NULL}, a default Daniell kernel with 
widths (9,9,9) is used.}

\item{segments}{Numeric value indicating the number of windowing segments to
use with \code{welchPSD} method.}

\item{...}{Additional arguments to be passed to spectral analysis functions, 
such as the \code{windowfun} option for \code{welchPSD}.}
}
\value{
List of wave parameters based on spectral methods.
\itemize{
  \item \code{h} Average water depth. Same units as input surface heights 
  (typically meters).

  \item \code{Hm0} Significant wave height based on spectral moment 0. Same
  units as input surface heights (typically meters).
  This is approximately equal to the average of the highest 1/3 of the waves.

  \item \code{Tp} Peak period, calculated as the frequency with maximum power 
  in the power spectrum. Units of seconds.

  \item \code{m0} Estimated variance of time series (moment 0).

  \item \code{T_0_1} Average period calculated as \eqn{m0/m1}, units seconds. Follows National 
  Data Buoy Center's method for average period (APD).

  \item \code{T_0_2} Average period calculated as \eqn{(m0/m2)^{0.5}}, units seconds. Follows 
  Scripps Institution of Oceanography's method for calculating average period 
  (APD) for their buoys.

  \item \code{EPS2} Spectral width parameter.

  \item \code{EPS4} Spectral width parameter.
}
}
\description{
Calculate ocean wave parameters using spectral analysis methods
}
\details{
Carries out spectral analysis of ocean wave height time series to estimate
common wave height statistics, including peak period, average period, 
and significant wave height.
}
\examples{
data(wavedata)
waveStatsSP(wavedata$swDepth.m, Fs = 4, method = 'spec.pgram', plot = TRUE)
}
\references{
Original MATLAB function by Urs Neumeier:  
http://neumeier.perso.ch/matlab/waves.html, based on code developed by Travis 
Mason, Magali Lecouturier and Urs Neumeier.
}
\seealso{
\code{\link{waveStatsZC}} for wave statistics determined using a 
zero-crossing algorithm.
}
