\name{GENEAread-package}
\alias{GENEAread-package}
\alias{GENEAread}
\docType{package}
\title{
GENEAread: a package to process binary accelerometer output files.
}
\description{
This is a package to process binary output files from the GENEA accelerometer data.  The main functions are:\cr\cr

read.bin\cr
stft\cr
epoch.apply\cr
}
\section{Main tasks performed}{
The main tasks performed by the package are listed below. The relevant topic contains documentation and examples for each.
\describe{
\item{Extraction of file header material}{is accomplished by \code{\link{header.info}}.}
\item{Input and downsampling of data}{is accomplished by \code{\link{read.bin}}.
}
\item{Selection of time intervals}{is accomplished via \code{\link{get.intervals}}.}
\item{Computation of epochal summaries}{is accomplished by \code{\link{epoch.apply}} and other functions documented therein.}
\item{Computation of STFT analyses}{is accomplished by \code{\link{stft}}.}
}
}
\section{Classes implemented}{
The package provides definitions and methods for the following S3 classes:
\describe{
\item{GRtime:}{Provides numeric storage and streamlined plotting for times. \code{\link{GRtime}}}
\item{AccData:}{Stores GENEA accelerometer data, allowing plotting, subsetting and other computation. \code{\link{AccData}}}
\item{VirtAccData:}{A virtual AccData object, for just-in-time data access via \code{\link{get.intervals}}.}
\item{stft:}{Processed STFT outputs, for plotting via \code{\link{plot.stft}}.}
}
}

\details{
\tabular{ll}{
Package: \tab GENEAread\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 19/06/2012\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Zhou Fang <zhou@activinsights.co.uk>
ActivInsights Ltd. <joss.langford@activinsights.co.uk>
}
\keyword{ package }
