\name{NEWS}
\title{News for Package \pkg{GENEAread}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 1.0}{
\itemize{
	\item Many changes and alterations and documentation.
	\item Renamed to GENEAread
	\item Version submitted to CRAN
}

\section{Changes in ReadGenea 2.5}{
\itemize{
	\item BUGFIX: Software should deal with multi-line comments correctly now.
	\item Package is now byte-compiled on installation - should lead to perfmance enhancement.
}

\section{Changes in ReadGenea 2.4}{
 \itemize{
	\item Implemented mmap-based reading of data - reading with mmap enabled should be much faster - note glitch exists with data recorded at very cold temperatures on existing firmware.
}

\section{Changes in ReadGenea 2.3}{
 \itemize{
	\item BUGFIX: Changed handing of timezones, times and dates should be more reasonable now
	\item Allowed 'day time' references to periods
	\item Implemented 'read from file' for virtual reads
}
}

\section{Changes in ReadGenea 2.2}{
 \itemize{
	\item Implemented some class features for matrix-like addressing of read databases
}
}

\section{Changes in ReadGenea 2.1}{
 \itemize{
    \item Added downsampling and block based loading of data (for very large datasets)    
    \item Added specification of start and end by proportions of data
    \item Some speed improvements/optimisations
    \item Experimental implementations of stft, GFLasso and FLasso.    
}
}

\section{Changes in ReadGenea 2.0-3}{
 \itemize{
    \item BUGFIX: Issue with reading of dates
}
}

\section{Changes in ReadGenea 2.0-2}{
  \itemize{
    \item Speed optimisations for read.bin/convert.hexstream
    \item Removed row labels on created object to reduce memory overheads
    \item NEWS file created
  }
}
