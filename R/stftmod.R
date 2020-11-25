
#' @name stft
#'
#' @title Computes Short Time Fourier Transforms
#'
#' @description Processes a dataset, creating an object contained processed time-frequency analyses. These can then be plotted.
#'
#' @param X The dataset to be processed.
#' @param start,end,length,time.format A specification for the segment to process, as in get.intervals.
#' @param type The type of STFT to compute.
#' @param mv.indices For type = "mv" or type = "sum", the indices to process and the order to process them in.
#' @param date.col logical. Whether the first column should be ignored and treated as a timestamp. If unset, is automatically chosen.
#' @param reassign logical. If TRUE, compute the time-reassigned STFT. For type c("mv", "sum"), this is done with the first coordinate in mv.indices.
#' @param plot.it logical. Whether to plot the STFT immediately when processing is complete, using the default plot.stft options.
#' @param ... Additional optional arguments to control the STFT computation. These are:\itemize{
#'   \item win: Window size in seconds for STFT computation. Increased window size mean better frequency resolution, but poorer time resolution. Defaults to 10 seconds.
#'   \item inc: Increment between successive time steps for processing. Defaults to win/2.
#'   \item coef: Number of fourier frequencies to compute. Small values will remove the higher frequencies from the processed object. Defaults to the maximum, win/2.
#'   \item wtype: String giving the name of a window function, providing coefficients for filtering before processing. "hanning.window" is the default, with "uniform.window" also available.
#'   \item freq: Sampling frequency of data set. If not given, is taken from X itself, or assumed to be 1 if unavailable.
#'   \item centre: If TRUE (Default), centre the data in each window before processing is done. Useful for avoiding excessively large DC offset coefficients in results.
#'   \item calc.null: If TRUE (Defaults to FALSE), compute a 'null' STFT by resampling the data completely, then doing a STFT.
#'   \item pvalues: If TRUE (Defaults to FALSE) Compute bootstrapped pvalues for each position by resampling within each window and applying a wilcox test.
#'   \item time: Allows the user to set an overriding timestamp vector to be used for processing.
#'   \item quiet: If TRUE, suppress output.
#' }
#'
#' @details This function accepts input in a variety of forms and computes short time fourier transforms to
#' extract frequency structure from the data.X may be an array, a vector, or an AccData object.
#' If date.col is TRUE, the first column of an array X would be used to determine timestamps.
#' Otherwise indices would be used. If date.col is not set, the function will attempt to determine whether
#' the first column is timestamp-like. The timestamp column is removed from X (and so not included in
#' consideration of mv.indices, for instance).
#' With vectors, the basic method is to compute the STFT by creating windows of size win seconds every
#' inc seconds, and computing the fourier transform. With multi-dimensional data and AccData,
#' processing is done on the dimensions that are in mv.indices, or the first three non-date columns
#' if that is unavailable. Three methods are possible:\itemize{
#'  \item 1. type = "mv": The one dimensional method is first applied to each of the chosen column indices. These are then collated by choosing, for each time-frequency combination, the maximum such value across each of the indices.
#'  \item 2. type = "svm": The SVM is computed first for each time step by computing the square rooted sum of squares. This is then dealt with using the one dimensional method.
#'  \item 3. type = "sum": As in "mv", the 1d method is applied. The square of the modulus of the result is then summed and square rooted.
#' }
#' If reassign is set, the time reassigned stft is also computed for the first element of mv.indices or the svm as appropriate, by using finite differencing. This gives potentially better resolution results for data with a clear signal component.
#'
#' @return A "stft" class object - a list with the following components:\itemize{
#'  \item call: The function call.
#'  \item type: Type of STFT computed.
#'  \item values: Mod of FFT computed, with each row corresponding to a specific time increment.
#'  \item increment,windowsize,centre,sampling.frequency: Various control parameters used in the computation.
#'  \item null.logmean,null.logsd: Log of the square rooted mean and standard deviation of the Mod FFT squared for the randomised data, if calc.null = TRUE.
#'  \item p.values: Wilcoxian pvalues, if pvalues = TRUE.
#'  \item principals: Principal frequencies.
#'  \item frequency: Frequencies at which FFT is computed.
#'  \item time: Timestamps for FFT windows.
#'  \item LGD: Local group delay matrix for reassigned STFT.
#'  \item CIF: Channelized instantaneous frequency matrix for reassigned STFT.
#'}
#'
#' @references Fulop, S.A. & Fitz, K. (2006). Algorithms for computing the time-corrected instantaneous frequency (reassigned) spectrogram, with applications J Acoustical Society of America 119(1), 360--371.
#' Nelson. D.J. (2001). Cross-spectral methods for processing speech J Acoustical Society of America 110(1), 2575-2592.
#'
#' @usage stft(X, start=0, end=1, length=NULL,  time.format = c("auto"),
#'             type = c("mv", "svm", "sum"), mv.indices,
#'             date.col,  reassign = TRUE, plot.it = FALSE,...)
#' @export
#'
#' @examples
#' \dontrun{
#'#Some artificial data
#'time = 1:5000
#'#sum of two sine curves at 0.3 Hz and 0.05 Hz
#'f1 = 0.3; f2 = 0.05
#'sin1 = sin(time * f1 * 2*pi)
#'sin2 = sin(time * f2 * 2*pi)
#'#add a bit of noise
#'signal = sin1 + sin2 + 1*rnorm(5000)
#'#non-reassigned
#'stft(signal, plot = TRUE, reassign = FALSE, win = 100)
#'#reassigned
#'stft(signal, plot = TRUE, reassign = TRUE, win = 100)
#'
#'#add a third component: varying frequency.
#'stft(signal + sin(cumsum(seq(f2, f1, length = 5000))*2*pi),
#'                  plot = TRUE, reassign = TRUE, win = 100)
#'
#'# Real data
#'binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]
#'
#'# Read in the entire file, calibrated
#'procfile<-read.bin(binfile)
#'# Default is mv
#'stft(procfile, plot.it = TRUE)
#'# Try sum?
#'stft(procfile, plot.it = TRUE, type = "sum", reassign = FALSE)
#'
#'# Just look at the last 50% of the data
#'stft(procfile, start = 0.5, plot.it = TRUE)
#'
#'# not reassigned, svm
#'stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE)
#'# a narrower 5 second window means better time resolution
#'stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, win = 5)
#'# choose increments so as not to overlap
#'stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, win = 5, inc = 5)
#'# uniform windows
#'stft(procfile, type = "svm", reassign = FALSE, plot.it = TRUE, wtype = "uniform.window")
#'# Svm, reassigned, quietly
#'obj = stft(procfile, type = "svm", quiet = TRUE)
#'plot(obj, cex = 3, showmax = FALSE, mode = "pval")
#'
#' #example code
#' plot(stft(subs(mag, 0.94,0.96), win = 1024, plot = F, coef = 512), zlog = T, log="y")
#' plot(stft(subs(mag, 0.7,8), win = 1024, plot = F, coef = 512), zlog = T, log="y")
#' plot(stft(subs(mag, 0.0001,0.005), win = 1024, plot = F, coef = 512), zlog = T)
#' plot(stft(subs(mag, 0.7,0.8), win = 1024, plot = F), zlog = T, log = "y")
#'
#' plot(stft(rep(1, 1000) +
#'       c(sin(1:500/ 10 * 2*pi), rep(0, 500)) +
#'       c(rep(0, 300),sin(1:500/ 20 * 2*pi), rep(0, 200)),
#'      freq = 1, plot.it = F), log="x")
#'
#' stft(sin(1:1000 / (1 +sqrt(1000:1)) * 2 * pi), freq = 1)
#' stft(rep(1, 1000) + sin(1:1000/ 10 * 2*pi), freq = 1)
#' }

stft <- function(X,
                 start = 0,
                 end = 1,
                 length = NULL,
                 time.format = c("auto"),
                 type = c("mv", "svm", "sum"),
                 mv.indices,
                 date.col,
                 reassign = TRUE,
                 plot.it = FALSE,...){
  type = match.arg(type)
  call <- match.call()
  if (is.list(X)){
  #if (!is.null(X$freq)) freq = X$freq

    X = get.intervals(X,
                      start,
                      end,
                      length,
                      time.format,
                      incl.date=TRUE,
                      simplify  = TRUE)
    if (missing(date.col)){
      date.col = TRUE
      }
  }

  if (length(dim(X)) < 2){
    X = matrix(X, ncol = 1)
    }

  #is first col date-like?
  if (missing(date.col)){

    if ((ncol(X)> 1) && ( X[1,1] >= 365 * 60*60*24)){
    print("Assuming first column is time.")
    date.col = T
    } else {
    date.col = F
    }
  }

  if (missing(mv.indices)){
    mv.indices <- 1:(min(ncol(X) - date.col, 3))
    }

  if (type == "svm"){
    if (date.col){
      obj1 = stftcalc(cbind(X[,1], sqrt(rowSums(X[,mv.indices + 1, drop = F]^2))), reassign = reassign,...)
    } else {
      obj1 = stftcalc( sqrt(rowSums(X[,mv.indices, drop = F]^2)),reassign = reassign, ...)
    }
  obj1$type = "svm"
  } else {
    ind = mv.indices[1]
    if (date.col) {ind = c(1,ind +1)}

    obj1 = stftcalc(X[,ind], reassign = reassign,...)

    if (length(mv.indices) > 1){

      if (type == "mv"){
        for (ind in mv.indices[-1] ) {

          if (date.col){ind = c(1,ind +1)}
          obj = stftcalc(X[,ind], reassign = F, ...)
          obj1$values = pmax(obj1$values, obj$values)
        }
      obj1$type = c("mv", paste(mv.indices, collapse=""))
      } else {
        for (ind in mv.indices[-1] ) {

          if (date.col){
            ind = c(1,ind +1)
            }
          obj = stftcalc(X[,ind], reassign = F, ...)
          obj1$values = sqrt(obj1$values^2 + obj$values^2)
        }
        obj1$type = c("sum", paste(mv.indices, collapse=""))
      }
    }
  }

  obj1$call = call
  if (plot.it){plot(obj1)}
  obj1
}

#' @title plot stft
#' @name plot.stft
#'
#' @description Processes a dataset, creating an object contained processed time-frequency analyses. These can then be plotted.
#'
#' @param x "stft" class object to be processed.
#' @param mode What should be plotted? \itemize{
#'     \item "decibels": log10 of FFT modulus
#'     \item "modulus": Raw FFT modulus
#'     \item "pvalue": P-value of each frequence's modulus assuming that window was in fact white noise of equal equal standard deviation
#' }
#' @param log For \code{log = "y"}, use a log scale on the y axis.
#' @param showmax Vector or logical. Compute and plot the principle frequency components?
#' @param median logical. If TRUE, smooth the STFT plot in the time direction with a running median.
#' @param xaxis logical. If TRUE, plot pretty time axes.
#' @param topthresh For finite values, crop plot for frequencies higher than this value, and show a summary plot up top.
#' @param reassign logical. Plot reassigned stft, if available?
#' @param xlim,ylim Parameters controlling axes limits of plot.
#' @param new logical. If TRUE, make a new plot. Otherwise overlay on to existing plot.
#' @param zlim.raw Raw values at which to threshold values for computation of heatmap colours.
#' @param zlim.quantile Quantile values at which to threshold values for computation of heatmap colours.
#' @param cex Size of points for reassigned STFT plotting.
#' @param col Vector of colours to be used for plotting.
#' @param ... Additional arguments to be passed to methods.
#'
#' @details STFT objects are created by the \code{\link{stft}} function. These methods print some useful summary statistics about them, and produce plots.
#' \code{mode} determines the type of plot. "decibel" and "modulus" work with the raw values, while "pvalue" conducts some degree of normalisation in each time window and so is perhaps more useful for data showing a large variation in sd across different points in time. If the \code{null.calc} was set in the original stft argument, that is used - otherwise, an Exponential distribution is fit to each window, and the pvalues computed from that.
#'
#' By default, the function uses some empirical quantile based colour thresholds designed to give somewhat reasonable and informative plots. This can be overridden, however, by setting different \code{zlim.raw} or \code{zlim.quantile} results. This can be useful for comparing two different datasets.
#'
#' Reassigned stft plots are constructed, by default, when they are available, and when the original was not a "mv" stft. Unlike the heatmap used in the usual stft plot, a 2d scatterplot is used instead. This means that if there are few data points, it can be advantageous to set a higher \code{cex} value for larger points and better display.
#'
#' With Accelerometer data, often the frequencies of interest are concentrated at the lower frequencies. Topthresh crops the frequency display to show only those frequencies. A summary plot is show on the top, to compensate. Choosing a grid of frequencies, this plot draws one line to represent the energies present in the signal at that particular frequency, and higher. Black lines are drawn for frequencies less than 2/3 the \code{topthresh}, red lines for 2/3 - 1 times \code{topthresh}, and blue lines for frequencies higher than \code{topthresh}. Alternative, set \code{log = "y"} to put frequencies on a log scale.
#'
#' @return These functions are run for their side effects.
#'
#' @importFrom grDevices col2rgb gray palette rgb
#' @importFrom graphics layout plot lines abline image points
#' @importFrom stats fft na.omit pexp runmed wilcox.test
#'
#' @seealso
#' \code{\link{stft}}, \code{\link{image.default}}
#'
#' @examples
#'  \dontrun{# Real data
#'  binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]
#'
#' #Read in the entire file, calibrated
#'  procfile<-read.bin(binfile)
#'  #Create stft object
#'  obj = stft(procfile, type = "svm", quiet = TRUE)
#'  #Look at it
#'  print(obj)
#'
#'  plot(obj, cex = 5)
#'  plot(obj, showmax = FALSE, cex = 5) #suppress principals
#'
#'  #pval plot
#'  plot(obj, mode = "pval", cex = 5)
#'  #disable reassigned stft
#'  plot(obj, mode = "pval", reassign = FALSE)
#'  #median smoothing
#'  plot(obj, mode = "pval", reassign = FALSE, median = TRUE)
#'  #log scale frequency, no top bar
#'  dev.new(); plot(obj, mode = "pval", reassign = FALSE, topthresh = Inf, log = "y")
#'}
#' @export

plot.stft <- function (x,
                       mode = c("decibels", "modulus", "pval"),
                       log = "",
                       showmax = TRUE,
                       median = FALSE,
                       xaxis = TRUE,
                       topthresh,
                       reassign = (!(is.null(x$LGD)) && !("mv" %in% x$type)),
                       ylim,
                       xlim,
                       new = TRUE,
                       zlim.raw,
                       zlim.quantile,
                       cex,
                       col = gray(63:0/63),
                       ...)
{
  xv <- x$values
  if (missing(cex)){
    cex = 5
    if (nrow(xv) > 500){
      cex = 2
    }
  }

  if (missing(topthresh)){
    topthresh = Inf
    if (x$sampling.frequency > 30){
      topthresh = 15
    }
  }

  if (median){
    xv = apply(xv,2, function(t) (runmed(t, k = 1 + 2 * min((length(t)-1)%/% 2, ceiling(0.005*length(t))) )))
  }

  mode = match.arg(mode)
  if (mode == "decibels"){
    xv = log(xv)
    if (missing(zlim.raw)){

      if (missing(zlim.quantile)){
        zlim.raw = c(median(xv), Inf)
        if (!is.null(x$null.logmean)){
          zlim.raw = c(x$null.logmean, Inf)
        }
      }
      else {
        zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
      }
    }

    if (length(zlim.raw) == 1){
      zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
    }
    xv = constrain(xv, zlim.raw[1], zlim.raw[2])
  }
  else if (mode == "pval"){
    xv = t(apply(xv, 1, function(t)  -log10(1-pexp(t^2, 1/mean(t^2)) )))
    if (missing(zlim.raw)){
      if (missing(zlim.quantile)){
        zlim.raw = c(0, 15)
      }
      else {
        zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
      }
    }
    if (length(zlim.raw) == 1){
      zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
    }
    xv = constrain(xv, zlim.raw[1], zlim.raw[2])
  }
  else {
    if (missing(zlim.raw)){
      if (missing(zlim.quantile)){
        zlim.raw = c(0, Inf)
      }
      else {
        zlim.raw = c(quantile(xv, zlim.quantile[1]), quantile(xv, zlim.quantile[2]))
      }
    }
    if (length(zlim.raw) == 1){
      zlim.raw = c(zlim.raw, zlim.raw + abs(zlim.raw)*0.0001)
    }
    xv = constrain(xv, zlim.raw[1], zlim.raw[2])
  }

  timegrid = x$times
  if (missing(ylim)){
    ylim = range( x$frequency)
  }
  if (missing(xlim)){
    xlim = range(timegrid)
  }
  frequency= x$frequency

  if (topthresh < Inf){
    if (new){
      layout(matrix(c(1,2,2,2), ncol = 1))
      par(mar = c(0,1,0,0))
      par(oma = c(5, 4, 4, 2) + 0.1)

      binwidth = ceiling(length(timegrid) /100)
      topind = 1:(floor(length(timegrid) / binwidth) * binwidth)

      #res = epoch.apply(x$values, epoch.size = binwidth, function(t) apply(t,2,median)) #apply(x$values, 2, function(t)  apply(matrix(t, binwidth), 2, median))
      res = apply(x$values[topind,], 2, function(t)  apply(matrix(t, binwidth), 2, median))
      timegridtop = matrix(timegrid[topind], binwidth)[1,]
      ind = ceiling(ncol(xv) * 1:20/20)
      ylim =  c(0, quantile(sqrt(rowSums(res^2)), 0.95)*1.1)

      if (xaxis){
        plot(  convert.time(timegridtop), sqrt(rowSums(res^2)) , type="l", xaxt = "n", ylim =ylim, ...)
        axis( 1, convert.time(timegridtop), labels = F)
        timegridtop = convert.time(timegridtop)
      }
      else {
        plot(  (timegridtop), sqrt(rowSums(res^2))  , type="l", xaxt="n", ylim = ylim,...)
        axis(1, pretty(timegridtop), labels = F)
      }
      for (k in ind){
        if (frequency[k] <= (topthresh *2/3)){
          colour = "black"
        }
        else if (frequency[k] %bt% c(topthresh * 2/3, topthresh)){
          colour = "red"
        }
        else {
          colour = "blue"
        }
        lines(timegridtop,sqrt(apply((res)[, k:ncol(res), drop=F]^2, 1, sum)), col=  colour)
        abline(h = 0)
      }
    }

    ylim[2] = min(ylim[2], topthresh)
    #xv[, which(frequency > topthresh)]
    xv = xv[,which(frequency <= topthresh)]
    #do top thresholding
  }

  if(log == "y"){
    frequency[1] = frequency[2]^2/frequency[3]
    frequency = c(frequency, tail(frequency,1)^2/tail(frequency,2)[1])
    ylim[1] = max(min(frequency), ylim[1])
  }
  if (reassign){
    frequency =  as.vector(x$CIF[,1:ncol(xv)])
    if (is.numeric(col)){
      colours = col2rgb(palette()[col] )
      colours = rgb(colours[1], colours[2], colours[3], alpha = 255*(conv01(as.vector(xv)))  , maxColorValue = 255)
    }
    else {
      colours = col[1+ (length(col) - 1) * ( conv01(as.vector(xv)))]
    }
  }
  if (xaxis){
    if (reassign){
      time = convert.time(rep(x$times, ncol(xv) )+ as.vector(x$LGD[,1:ncol(xv)] ))
      if (new){
        plot( time, frequency,pch= ".", cex = cex , col = colours , log = log,ylim = ylim , xlim = convert.time(xlim), ...)
      }
      else{
        points ( time, frequency,pch= ".", cex = cex , col = colours , ylim = ylim ,  ...)
      }
      #####
    }
    else {
      time = timegrid
      plot(convert.time(  seq(min(time), max(time), len = 20) ), rep(1,20), col=0, xlab = "time", ylab = "frequency", ylim = ylim, xlim = convert.time(xlim), xpd = NA, log =log)
      #	par(new = T)
      image( x = convert.time(time) , y = frequency[1:ncol(xv) ],   z=xv, col=col, log = log, xaxt = "n", ylim = ylim,xlim = convert.time(xlim), add= T,...)
    }
  }
  else {
    if (reassign){
      time = (rep(x$times, ncol(xv) )+ as.vector(x$LGD[,1:ncol(xv)] ))
      if (new){
        plot( time, frequency,pch= ".", cex = cex , col = colours , log = log,ylim = ylim , xlim = xlim, ...)
      }
      else {
        points( time, frequency,pch= ".", cex = cex , col = colours , ylim = ylim , ...)
      }
      #######
    }
    else {
      time = timegrid
      image( x = time , y = frequency[1:ncol(xv)],   z=xv, col=col, log = log, ylim = ylim,xlim=xlim,...)
    }
  }
  axis(2, pretty(constrain(frequency, ylim[1], ylim[2]), min(floor(topthresh), floor(max(frequency))) ), labels = F, tcl = -0.2)
  if (as.numeric(showmax) > 0){
    #points ( time, x$principals, col=2 * (rowMeans(xv) > 1 * x$null.logmean)  , pch=".", cex = 3)
    pseudolines(timegrid, x$principals, col = 2, pch = ".", lwd = 2, cex = 2)
  }
  if (as.numeric(showmax) > 1){
    pseudolines(timegrid, frequency[ apply(x$values, 1, function(t) which.max(replace(t, which.max(t), -Inf)))], col=3, pch = ".", lwd = 2, cex = 2)
  }
}

#' @name stftcalc
#'
#' @title calculate the fft components
#'
#' @description  stftcalc calculating the short time fourier transform component
#'
#' @param x The dataset to be processed
#' @param win Window size in number of recordings to use.
#' @param inc increment size from one window to the next
#' @param coef Number of fourier frequencies to compute. Small values will remove the higher frequencies from the processed object. Defaults to the maximum, win/2.
#' @param wtype Window type for the STFT calculation.
#' @param freq frequency of dataset.
#' If missing and length(time) > 1, freq =(length(time) -1) /( max(time) - min(time) )
#' If just missing, freq = 1
#' @param centre If TRUE (Default), centre the data in each window before processing is done.
#' Useful for avoiding excessively large DC offset coefficients in results.
#' @param calc.null If TRUE (Defaults to FALSE), compute a 'null' STFT by resampling the data completely, then doing a STFT.
#' @param pvalues If TRUE (Defaults to FALSE) Compute bootstrapped pvalues for each position by resampling within each window and applying a wilcox test.
#' @param time Allows the user to set an overriding timestamp vector to be used for processing.
#' @param reassign logical. If TRUE, compute the time-reassigned STFT. For type \%in\% c("mv", "sum"), this is done with the first coordinate in mv.indices.
#' @param quiet If TRUE, suppress output.
#'
#' @importFrom utils object.size setTxtProgressBar txtProgressBar tail
#'
#' @keywords internal

stftcalc <- function(X,
                     win   = 10,
                     inc   =  win/2,
                     coef  = Inf,
		                 wtype = "hanning.window",
		                 freq ,
		                 centre = T,
		                 calc.null = F ,
		                 pvalues = F,
		                 time = NULL,
		                 reassign = T ,
		                 quiet = F){

  call = match.call()
  if (length(dim(X)) ==2) {
    if (is.null(time)){
      time = X[,1]
    }
    X = X[,2]
  }

  if ((length(time) >1) && (missing(freq))){
    freq =(length(time) -1) /( max(time) - min(time) )
  }
  if (missing(freq)){
    freq = 1
  }
  # Initial variables
  inc0 = inc
  win0 = win
  coef0 = coef
  win = round(win * freq)
  inc = max(round(inc * freq),1)
  coef = min(coef, floor(win/2))

  Xdel = shift(X, c(1,0), fill = "edge")
  numcoef <- 2*coef

  if (win < numcoef){
  	win <- numcoef
    if (!quiet){
      cat ("stft: window size adjusted to", win, ".\n")
    }
  }
  numwin <- trunc ((length(X) - win) / inc)

  ## compute the windows coefficients
  wincoef <- eval(parse(text = wtype))(win)

  ## create a matrix Z whose columns contain the windowed time-slices
  pval = rep(0, numwin+1)
  z <- matrix (0, numwin + 1, win)
  y <- matrix(0, numwin+1, win)
  ydel <- matrix(0, numwin+1, win)
  st <- 1

  if (!quiet){
    pb <- txtProgressBar(min = 0, max = 100,style=1)
  }

  for (i in 0:numwin){
  	z[i+1, 1:win] <- (X[st:(st+win-1)] - mean(X[st:(st+win - 1)])* centre) * wincoef
  	y[i+1,] <- fft(z[i+1,] )

    if (reassign){
    	z[i+1, 1:win] <- (Xdel[st:(st+win-1)] - mean(Xdel[st:(st+win - 1)])* centre) * wincoef
    	ydel[i+1,] <- fft(z[i+1,] )
    }

    if (pvalues){
      temp = sample(X[st:(st + win - 1)])
      temp = (temp - mean(temp) * centre)*wincoef
      temp = Mod(fft(temp))[1:coef]^2
      pval[i+1] = wilcox.test( (Mod(y[ i+1,])^2 - mean(Mod(y[ i+1,])^2))^2, (temp - mean(temp))^2)$p.value
    }
  	st <- st + inc
    if (!quiet){
      setTxtProgressBar(pb, 90* i/(numwin ))
    }
  }

  if (reassign){
    yfreqdel = cbind(y[, win],y[, 2: win - 1])# t(apply(y, 1, function(t) shift(t, 1, fill = "loop")))
    ydel = Arg(y*Conj(ydel)) *(freq/(2*pi))
    yfreqdel = -(win/(2*pi*freq)) * Arg(y * Conj(yfreqdel)) + win/(2*freq)
  } else {
    yfreqdel = NULL
  }

  null.logmean = NULL
  null.logsd = NULL

  if (calc.null){
    tmpdat = stft(sample(X),
                  win = win0,
                  inc= inc0,
                  coef=coef0,
                  wtype=wtype,
                  freq = freq,
                  centre = T,
                  calc.null = F,
                  quiet= T)
    null.logmean = log(sqrt(mean((tmpdat$values)^2)))
    null.logsd = log(sd(tmpdat$values))
  }

  if (!quiet){
    setTxtProgressBar(pb, 100)
  }

  if (is.null(time)){
      Y <- list (values = cbind(Mod(y[,1]),
                          2*Mod(y[,(2):coef])),
                          windowsize = win,
                          increment = inc,
                          windowtype = wtype,
                          centre = centre,
                          sampling.frequency = freq,
                          null.logmean = null.logmean,
                          null.logsd = null.logsd,
                          principals = (freq * (1:coef  - 1 ) / win)[apply( Mod(y[,(1):coef]),1, which.max)],
                          frequency = (freq * (1:coef  - 1 ) / win),
                          times =  (win/2 +  inc * 0:(nrow(y) - 1))/(freq),
                          p.values = pval,
                          LGD = yfreqdel[,1:coef],
                          CIF = ydel[,1:coef] )
  } else {
    if (length(time) == 1){
      times = (time  + (win/2 +  inc * 0:(nrow(y) - 1))/freq)
    } else {
      times = time[win/2 + inc* 0:(nrow(y) - 1) ]
    }

    Y <- list(values = cbind(Mod(y[,1]),
                             2*Mod(y[,(2):coef])),
                             windowsize = win,
                             increment = inc,
                             windowtype = wtype,
                             centre = centre,
                             sampling.frequency = freq,
                             null.logmean = null.logmean,
                             null.logsd = null.logsd,
                             principals = (freq * (1:coef  - 1 ) / win)[apply( Mod(y[,(1):coef]),1, which.max)],
                             frequency = (freq * (1:coef  - 1 ) / win),
                             times = times,
                             p.values = pval,
                             LGD = yfreqdel[,1:coef],
                             CIF = ydel[,1:coef]  )
  }

  if (!quiet){
    close(pb)
  }
  Y$call = call
  class(Y) <- "stft"
  return(Y)
  }

#' @title Hanning Window
#'
#' @description A hanning window used by the STFT function
#'
#' @param n number of points inside the window
#'
#' @keywords Internal

hanning.window = function (n) {
  if (n == 1)
      c <- 1
  else {
      n <- n - 1
      c <- 0.5 - 0.5 * cos(2 * pi * (0:n)/n)
  }
  return(c)
}

#' @title Uniform Window
#'
#' @description A uniform window used by the STFT function
#'
#' @param n number of points inside the window
#'
#' @keywords Internal

uniform.window = function(n){
  rep(1, n)
}

#' @title subs
#'
#' @description subsets a proportion of the dataset, or a certain length of the dataset starting at a
#' specific proportion position
#'
#' @param x object to subset
#' @param a Start point to separate the object
#' @param b End point to separate the object
#'
#' @keywords internal

subs <- function(x, a,b){
  len = length(x)
  if (a > 1){
  return(x[a : (a+b - 1)])
  } else if (b > 1) {
  return( x [ floor(a * len): (floor(a*len) + b - 1)])
  } else {
  return (x[floor(a*len) : floor(b*len)])
  }
}

#' @title getfreqs
#'
#' @description gets fft components corresponding to frequencies
#'
#' @param x stft object
#' @param frequencies frequencies range
#'
#' @importFrom stats fft
#'
#' @keywords internal

getfreqs = function(x, frequencies){
  n = length(x)
  fftobj = fft(x)
  frequencies = c(frequencies, n - frequencies) + 1
  frequencies = frequencies[which(frequencies != n+1)]
  fftobj = replace(fftobj, (1:n)[ - frequencies], 0)
  return(Re(fft(fftobj, inverse=T))/n)
}

#' @title print.stft
#'
#' @description print the stft object
#'
#' @param x stft object
#' @param ... ignored
#'
#' @importFrom utils tail
#'
#' @export

print.stft = function(x,...){
  cat("STFT object:\n")
  cat(format.GRtime(x$times[1],  format = "%y-%m-%d %H:%M:%OS3 (%a)")," to ", format.GRtime((tail(x$times,1)), format = "%y-%m-%d %H:%M:%OS3 (%a)"), "\n")
  cat(nrow(x$values), "increments of" , round(x$increment/x$sampling.freq, 3), "s \n")
  cat("Window size: " , x$windowsize, "(", round(x$windowsize/x$sampling.frequency, 3), "s ) -> f resolution: ", round(x$frequency[2],3), "Hz\n")
  if ("svm" %in% x$type) cat("[SVM]")
  if ("mv" %in% x$type) cat("[MV-", x$type[2], "]")
  if ("sum" %in% x$type) cat("[SUM-", x$type[2], "]")
  if (!is.null(x$LGD)) cat ("[Reassign]")
  cat("\n------ \n")
  cat("{" ,format(x$call), "}\n")
}

