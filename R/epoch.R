#' @name epoch
#' @aliases epoch
#' @aliases apply.epoch
#' @aliases mean.epoch
#' @aliases sd.epoch
#' @aliases median.epoch
#' @aliases mad.epoch
#' @aliases autocor.epoch
#' @aliases quantile.epoch
#'
#' @title Compute epochal summary statistics.
#'
#' @description Computes epochal summary statistics for an "AccData" object, matrix, or vector, and collates into a matrix or vector.
#'
#' @param x The object to compute statistics for. Can be an "AccData" object, a matrix, or a vector.
#' @param epoch.size Numeric giving intervals to consider and aggregate. For "AccData" x taken as seconds. Otherwise, considered as rows, or as individual readings.
#' @param incl.date logical. If TRUE, include a column of times or original indices with the results.
#' @param FUN A function to be applied to each epoch.
#'
#' @details These functions compute epochal summary statistics for "AccData" objects, matrices and vectors.
#'
#' apply.epoch is the general function - according to the size of epoch.size, it splits up the x into collections
#' of consecutive rows, each with the same size. These are then successively supplied to FUN as its first argument.
#' If the result of FUN is a single value, then the results are concatenated into a vector output.
#' Otherwise, an array is formed with each row corresponding to a single epochal group. For AccData, the sampling
#' frequency of the dataset is used to interpret the epoch size in seconds. Otherwise, the raw record indices are used.
#' If incl.date is set, the original timestamp vector of the data, or the original indices, are downsampled and
#' included as the first column of the output.
#'
#' The remaining functions are wrappers that compute various commonly useful statistics -- in particular, applied to
#' "AccData" objects and arrays, they by default compute the epochal SVM mean, standard deviation, median, median
#' absolute deviation, and autocorrelation, and sample quantiles respectively. (Arrays are treated as each column
#' representing the x, y, and z components respectively.) Applied to vector input, processing will occur without the
#' SVM calculation. This behaviour may be overridden by the sqrt setting, which will force the function to use the
#' squared (default for arrays and "AccData") or original unit (default for vectors) values in the statistical analysis.
#'
#' @return A vector or array giving the computed epochal summaries. With incl.date = TRUE,
#' the result is given as a data.frame suitable for plotting.
#'
#' @importFrom stats acf mad median sd
#'
#' @examples
#' \dontrun{
#' dat <- read.bin(system.file("binfile/TESTfile.bin", package = "GENEAread")[1]
#' , calibrate = TRUE)
#'
#' #look for the epochs that exceed a certain threshold 50% of the time
#' plot(apply.epoch( dat, epoch.size = 3 ,
#'                   FUN = function(t) mean(abs(svm(t) -1)>0.2)> 0.5 ), type = "l")
#'
#' plot(dat[,1], svm(dat), log = "y", pch = ".")
#' lines(mean.epoch(dat, incl.date = TRUE), lwd = 2)
#' lines(mean.epoch(dat, epoch.size = 30, incl.date = TRUE), col = 2, lwd = 2)
#' # This should give all the same results, but by a different way
#' lines(apply.epoch(dat, epoch.size = 30,
#'                   FUN = function(A) mean(svm(A, FALSE)), incl.date = TRUE), col = 3)
#' epsize = 30
#' lines(apply.epoch(dat, epoch.size = epsize,
#'                   FUN = function(t) median(t[,1])),
#'                   apply.epoch(dat, epoch.size = epsize,
#'                   FUN = function(A) mean(svm(A, FALSE))), col = 4)
#' #note this is different
#' lines(apply.epoch(dat, epoch.size = epsize,
#'                   FUN = function(t) median(t[,1])),
#'                   apply.epoch(dat, epoch.size = epsize,
#'                               FUN = function(A) mean(svm(A, sqrt = TRUE)))^2,
#'                               col = 5)
#'
#' #plot some statistics
#' par(mfrow = c(5,1), mar = c(1,4.5,1,1))
#' plot(sd.epoch(dat), type="l")
#' plot(median.epoch(dat), type= "l")
#' plot(mad.epoch(dat), type= "l")
#' plot(autocor.epoch(dat), type= "l")
#' tmp = quantile.epoch(dat, quantiles= c(0.1, 0.25, 0.5, 0.75, 0.9)); matplot(tmp, type = "l")
#'}
#' @export

apply.epoch <- function(x, epoch.size=10, incl.date = FALSE, FUN){
  sampling.freq = 1
  if (length(dim(x)) < 2){x = matrix(x, ncol = 1)}
  ind = 1:nrow(x)
  if (inherits(x, "AccData")){
    sampling.freq = x$freq
    times = x[,1]
    x = x$data.out#[,2:4]
  } else {
    times = ind
  }
  epoch.size = floor(epoch.size* sampling.freq)
  if (length(FUN(x[1:epoch.size,])) > 1){
    x = bapply(ind, epoch.size, function(t) FUN(x[t,]))
  } else {
    x = bapply.basic(ind, epoch.size, function(t) FUN(x[t,]))
  }
  if (incl.date){
    x = data.frame(time = times[seq(1, length(times)-epoch.size+1, by = epoch.size) + ceiling(epoch.size/2)],value = x)
  }
  x
}

#' @export

mean.epoch<- function(x, epoch.size=10,
                      incl.date = FALSE, sqrt,...) {

  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
    apply.epoch(x, epoch.size = epoch.size,
                FUN = function(x) mean(svm(x, sqrt)),
                incl.date= incl.date)
}

#' @method sd epoch
#'
#' @export

sd.epoch<- function(x, epoch.size=10, incl.date = FALSE, sqrt,...){
  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
  apply.epoch(x, epoch.size = epoch.size, FUN = function(x) sd(svm(x, sqrt)), incl.date= incl.date)
  }

#' @export

median.epoch<- function(x, na.rm = TRUE, epoch.size=10, incl.date = FALSE, sqrt,...){
  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
  apply.epoch(x, epoch.size = epoch.size, FUN = function(x) median(svm(x, sqrt), na.rm = na.rm), incl.date= incl.date)}


#' @method mad epoch
#'
#' @export

mad.epoch <- function(x, epoch.size=10, incl.date = FALSE, sqrt,...){
  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
  apply.epoch(x, epoch.size = epoch.size, FUN = function(x) mad(svm(x, sqrt)), incl.date= incl.date)}

#' @method acf epoch
#'
#' @export

acf.epoch<- function(x, epoch.size=10, lag = 1,
                         type = c("correlation", "covariance", "partial"),
                         incl.date = FALSE, sqrt,... ){
  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
  apply.epoch(x, epoch.size = epoch.size, FUN = function(x) acf(svm(x, sqrt), type = type , lag.max = max(lag), plot = F)$acf[lag +1], incl.date= incl.date)}



#' @method quantile epoch
#'
#' @export

quantile.epoch <- function(x, epoch.size = 10, quantiles= c(0.1, 0.25, 0.5, 0.75, 0.9),
                           incl.date = FALSE, sqrt,... ){
  if (missing(sqrt)) sqrt = (length(dim(x)) < 2)
  apply.epoch(x, epoch.size = epoch.size, FUN = function(x) quantile(svm(x, sqrt), prob = quantiles, name = F), incl.date= incl.date)}



bapply.basic <- function(X, k, FUN){
  res = rep(0, floor(length(X) / k))
  for (i in 1:floor(length(X)/k)){
    res[i] = FUN(X[ (i-1)*k + 1:k])
  }
  return(res)
  }

bapply <- function(X, k, FUN){
  dimout = length(FUN(X[1:k]))
  res = matrix(0, dimout, floor(length(X) / k))
  for (i in 1:floor(length(X)/k)){
    res[(i-1)* dimout + 1:dimout] = FUN(X[ (i-1)*k + 1:k])
  }
  return(t(res))
}



