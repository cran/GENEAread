\name{stft}
\alias{stft}
\title{
Computes Short Time Fourier Transforms
}
\description{
Processes a dataset, creating an object contained processed time-frequency analyses. These can then be plotted.
}
\usage{
stft(X, start=0, end=1, length=NULL,  time.format = c("auto"), 
    type = c("mv", "svm", "sum"), mv.indices, date.col, 
    reassign = TRUE,plot.it = FALSE,\dots)
}
\arguments{
  \item{X}{The dataset to be processed.}
  \item{start, end, length, time.format}{A specification for the segment to process, as in \code{\link{get.intervals}}.}
  \item{type}{The type of STFT to compute.}
  \item{mv.indices}{For \code{type = "mv"} or \code{type = "sum"}, the indices to process and the order to process them in.}
  \item{date.col}{logical. Whether the first column should be ignored and treated as a timestamp. If unset, is automatically chosen.}
  \item{reassign}{logical. If TRUE, compute the time-reassigned STFT. For \code{type \%in\% c("mv", "sum")}, this is done with the first coordinate in \code{mv.indices}.}
  \item{plot.it}{logical. Whether to plot the STFT immediately when processing is complete, using the default \code{plot.stft} options.}
  \item{\dots}{Additional optional arguments to control the STFT computation. These are: \cr\cr

  \code{win}: Window size in seconds for STFT computation. Increased window size mean better frequency resolution, but poorer time resolution. Defaults to 10 seconds.\cr

  \code{inc}: Increment between successive time steps for processing. Defaults to \code{win/2}. \cr

  \code{coef}: Number of fourier frequencies to compute. Small values will remove the higher frequencies from the processed object. Defaults to the maximum, \code{win/2}.\cr

  \code{wtype}: String giving the name of a window function, providing coefficients for filtering before processing. "hanning.window" is the default, with "uniform.window" also available.\cr

  \code{freq}: Sampling frequency of data set. If not given, is taken from X itself, or assumed to be 1 if unavailable.\cr

  \code{center}: If TRUE (Default), center the data in each window before processing is done. Useful for avoiding excessively large DC offset coefficients in results.\cr

  \code{calc.null}: If TRUE (Defaults to FALSE), compute a 'null' STFT by resampling the data completely, then doing a STFT.\cr

  \code{pvalues}: If TRUE (Defaults to FALSE) Compute bootstrapped pvalues for each position by resampling within each window and applying a wilcox test.\cr

  \code{time}: Allows the user to set an overriding timestamp vector to be used for processing.\cr

  \code{quiet}: If TRUE, suppress output.
}
}
\details{
This function accepts input in a variety of forms and computes short time fourier transforms to extract frequency structure from the data.

X may be an array, a vector, or an AccData object. If date.col is TRUE, the first column of an array X would be used to determine timestamps. Otherwise indices would be used. If date.col is not set, the function will attempt to determine whether the first column is timestamp-like. The timestamp column is removed from X (and so not included in consideration of \code{mv.indices}, for instance).

With vectors, the basic method is to compute the STFT by creating windows of size \code{win} seconds every \code{inc} seconds, and computing the fourier transform. With multi-dimensional data and AccData, processing is done on the dimensions that are in \code{mv.indices}, or the first three non-date columns if that is unavailable. Three methods are possible:

1. \code{type = "mv"}: The one dimensional method is first applied to each of the chosen column indices. These are then collated by choosing, for each time-frequency combination, the maximum such value across each of the indices.

2. \code{type = "svm"}: The SVM is computed first for each time step by computing the square rooted sum of squares. This is then dealt with using the one dimensional method.

3. \code{type = "sum"}: As in "mv", the 1d method is applied. The square of the modulus of the result is then summed and square rooted.

If \code{reassign} is set, the time reassigned stft is also computed for the first element of \code{mv.indices} or the svm as appropriate, by using finite differencing. This gives potentially better resolution results for data with a clear signal component.
}

\section{Acknowledgements}{
The initial implementation of this function is based on that in package \code{e1071}. The reassigned STFT implementation is due to Nelson (2001), via Fulop (2006).
}
\references{
Fulop, S.A. & Fitz, K. (2006). Algorithms for computing the time-corrected instantaneous frequency (reassigned) spectrogram, with applications \emph{J Acoustical Society of America} \bold{119(1)}, 360--371.

Nelson. D.J. (2001). Cross-spectral methods for processing speech \emph{J Acoustical Society of America} \bold{110(1)}, 2575-2592.
}

\value{
A "stft" class object - a list with the following components:
\item{call}{The function call.}
\item{type}{Type of STFT computed.}
\item{values}{Mod of FFT computed, with each row corresponding to a specific time increment.}
\item{increment, windowsize, center, sampling.frequency}{Various control parameters used in the computation.}
\item{null.logmean, null.logsd}{Log of the square rooted mean and standard deviation of the Mod FFT squared for the randomised data, if \code{calc.null = TRUE}.}
\item{p.values}{Wilcoxian pvalues, if \code{pvalues = TRUE}.}
\item{principals}{Principal frequencies.}
\item{frequency}{Frequencies at which FFT is computed.}
\item{time}{Timestamps for FFT windows.}
\item{LGD}{Local group delay matrix for reassigned STFT.}
\item{CIF}{Channelized instantaneous frequency matrix for reassigned STFT.}
}
\seealso{
\code{\link{plot.stft}}, \code{\link{get.intervals}}
}
\examples{
#Some artificial data
time = 1:5000
#sum of two sine curves at 0.3 Hz and 0.05 Hz
f1 = 0.3; f2 = 0.05
sin1 = sin(time * f1 * 2*pi)
sin2 = sin(time * f2 * 2*pi)
#add a bit of noise
signal = sin1 + sin2 + 1*rnorm(5000)
#non-reassigned
stft(signal, plot = TRUE, reassign = FALSE, win = 100)
#reassigned
stft(signal, plot = TRUE, reassign = TRUE, win = 100)

#add a third component: varying frequency.
stft(signal + sin( cumsum(seq(f2, f1, length = 5000))*2*pi), plot = TRUE, reassign = TRUE, win = 100)

#Real data
binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]

#Read in the entire file, calibrated
procfile<-read.bin(binfile)
#Default is mv
stft(procfile, plot.it = TRUE)
#Try sum?
stft(procfile, plot.it = TRUE, type = "sum", reassign = FALSE)

#Just look at the last 50% of the data
stft(procfile, start = 0.5, plot.it = TRUE)

#not reassigned, svm
stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE)
#a narrower 5 second window means better time resolution
stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, win = 5)
#choose increments so as not to overlap
stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, win = 5, inc = 5)
#uniform windows
stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, wtype = "uniform.window")

#Svm, reassigned, quietly
obj = stft(procfile, type = "svm", quiet = TRUE)
plot(obj, cex = 3, showmax = FALSE, mode = "pval")
}
\keyword{ts}
