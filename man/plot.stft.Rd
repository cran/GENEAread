\name{plot.stft}
\alias{plot.stft}
\alias{print.stft}
\title{
Plots and prints Short Time Fourier Transforms
}
\description{
Processes a dataset, creating an object contained processed time-frequency analyses. These can then be plotted.
}
\usage{
\method{plot}{stft}(x, mode = c("decibels", "modulus", "pval"), log = "", showmax = TRUE, median = FALSE, xaxis = TRUE, topthresh, reassign = (!(is.null(x$LGD)) && !("mv" \%in\% x$type)), ylim, xlim,new = TRUE, zlim.raw,zlim.quantile, cex, col = gray (63:0/63),\dots)
\method{print}{stft}(x, \dots)
}
\arguments{
  \item{x}{"stft" class object to be processed.}
  \item{mode}{What should be plotted?\cr\cr
  "decibels": log10 of FFT modulus\cr
  "modulus": Raw FFT modulus\cr
  "pvalue": P-value of each frequence's modulus assuming that window was in fact white noise of equal equal standard deviation}
  \item{log}{For \code{log = "y"}, use a log scale on the y axis.}
  \item{showmax}{Vector or logical. Compute and plot the principle frequency components?}
  \item{median}{logical. If TRUE, smooth the STFT plot in the time direction with a running median.}
  \item{xaxis}{logical. If TRUE, plot pretty time axes.}
  \item{topthresh}{For finite values, crop plot for frequencies higher than this value, and show a summary plot up top.}
  \item{reassign}{logical. Plot reassigned stft, if available?}
  \item{xlim, ylim}{Parameters controlling axes limits of plot.}
  \item{new}{logical. If TRUE, make a new plot. Otherwise overlay on to existing plot.}
  \item{zlim.raw}{Raw values at which to threshold values for computation of heatmap colours.}
  \item{zlim.quantile}{Quantile values at which to threshold values for computation of heatmap colours.}
  \item{cex}{Size of points for reassigned STFT plotting.}
  \item{col}{Vector of colours to be used for plotting.}
  \item{\dots}{Additional arguments to be passed to methods.}
}
\details{
STFT objects are created by the \code{\link{stft}} function. These methods print some useful summary statistics about them, and produce plots.

\code{mode} determines the type of plot. "decibel" and "modulus" work with the raw values, while "pvalue" conducts some degree of normalisation in each time window and so is perhaps more useful for data showing a large variation in sd across different points in time. If the \code{null.calc} was set in the original stft argument, that is used - otherwise, an Exponential distribution is fit to each window, and the pvalues computed from that.

By default, the function uses some empirical quantile based colour thresholds designed to give somewhat reasonable and informative plots. This can be overridden, however, by setting different \code{zlim.raw} or \code{zlim.quantile} results. This can be useful for comparing two different datasets.

Reassigned stft plots are constructed, by default, when they are available, and when the original was not a "mv" stft. Unlike the heatmap used in the usual stft plot, a 2d scatterplot is used instead. This means that if there are few data points, it can be advantageous to set a higher \code{cex} value for larger points and better display.

With Accelerometer data, often the frequencies of interest are concentrated at the lower frequencies. Topthresh crops the frequency display to show only those frequencies. A summary plot is show on the top, to compensate. Choosing a grid of frequencies, this plot draws one line to represent the energies present in the signal at that particular frequency, and higher. Black lines are drawn for frequencies less than 2/3 the \code{topthresh}, red lines for 2/3 - 1 times \code{topthresh}, and blue lines for frequencies higher than \code{topthresh}. Alternative, set \code{log = "y"} to put frequencies on a log scale.
}

\value{
These functions are run for their side effects.
}
\seealso{
\code{\link{stft}}, \code{\link{image.default}}
}
\examples{

#Real data
binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]

#Read in the entire file, calibrated
procfile<-read.bin(binfile)
#Create stft object
obj = stft(procfile, type = "svm", quiet = TRUE)
#Look at it
print(obj)

plot(obj, cex = 5)
plot(obj, showmax = FALSE, cex = 5) #suppress principals

#pval plot
plot(obj, mode = "pval", cex = 5)
#disable reassigned stft
plot(obj, mode = "pval", reassign = FALSE) 
#median smoothing
plot(obj, mode = "pval", reassign = FALSE, median = TRUE) 
#log scale frequency, no top bar
dev.new(); plot(obj, mode = "pval", reassign = FALSE, topthresh = Inf, log = "y") 
}
\keyword{hplot}

