% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detrendHeight.R
\name{detrendHeight}
\alias{detrendHeight}
\title{Remove trend from a time series}
\usage{
detrendHeight(pt)
}
\arguments{
\item{pt}{A vector of numeric values to be detrended}
}
\value{
A list containing the following:

\code{pt} A vector of detrended values 

\code{seg_len} The segment length used 

\code{h} The mean height of the water column

\code{trend} A two element vector of the intercept and slope
 from the linear regression.
}
\description{
Fits a straight line to a vector of values using lm(), and uses the
regression coefficients to subtract off the linear trend from the values.
}
\details{
Typically this is used to remove a tidal trend from a ocean 
surface height time series before attempting to calculate statistics for 
waves, which are variations above and below the mean surface height. 
Returns a series of residuals around the linear trend, in the original units
of measurement (typically meters).
}
\examples{
data(wavedata)
detrended <- detrendHeight(wavedata$swDepth.m)
pt <- detrended[['pt']]
plot(pt, type = 'l')
abline(h = 0)
}
