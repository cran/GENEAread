
#' @name Convert.time
#'
#' @title Convert.time
#'
#' @description extends time display from package chron to use h:m:s for > 1 day times.
#'
#' @param x Object to process. For convert.time, must be numeric.
#' @param format A character string indicating the form of output. See \code{\link[base]{strptime}} for details. If
#' NULL, will be automatically chosen.
#'
#' @details convert.time converts numerics to GRtime objects. The format argument allows a format string to be attached
#' specifying the default format to display in. as.GRtime is a wrapper to convert.time, that when supplied with
#' character input, coerces the value first to numeric using parse.time.
#'
#' @export

convert.time = function(x, format = NULL){
  if (!inherits(x, "GRtime")){class(x) = c("GRtime",class(x))}
  attr(x, "format") = format
  x
}
