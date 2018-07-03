#' @name GRtime
#' @aliases GRtime
#' @aliases c.GRtime
#' @aliases as.GRtime
#' @aliases format.GRtime
#' @aliases Axis.GRtime
#' @aliases axis.GRtime
#' @aliases Ops.GRtime
#' @aliases pretty.GRtime
#' @aliases print.GRtime
#'
#' @title Date time handling for the GENEAread package.
#'
#' @description Stores date time data as a numeric, with facility for pretty printing and axis commands.
#'
#' @param x Object to process. For \code{convert.time}, must be numeric. For \code{as.GRtime} may be
#' numeric or character. For \code{format.GRtime}, a GRtime object, or a numeric.
#' @param format A character string indicating the form of output. See \code{\link{strptime}} for details.
#' If NULL, will be automatically chosen.
#' @param ... Additional arguments to be passed to \code{\link{parse.time}}, \code{\link{as.numeric}},
#' \code{\link{format.POSIXct}}, \code{\link{axis}}, \code{\link{pretty.POSIXt}}.
#'
#' @details
#' The GRtime class handles dates and times for the GENEAread class. The class treats dates as numerics
#' denoting seconds since the UNIX epoch, with potentially a string attached specifying the format to print in.
#'  Unlike \code{POSIXct}, we avoid some of the processing, especially with respect to time zones, and allow
#'  some more flexibility in time computation and display. A range of operators are defined.
#'
#' \code{convert.time} converts numerics to GRtime objects. The \code{format} argument allows a format
#' string to be attached specifying the default format to display in. \code{as.GRtime} is a wrapper to
#' \code{convert.time}, that when supplied with character input, coerces the value first to numeric using
#' \code{parse.time}.
#'
#' \code{format.GRtime} formats GRtime objects for pretty printing. If \code{format} is provided as argument,
#'  that is used. Else, if the \code{format} attribute is set on \code{x}, that is used. Finally, if formats
#'  are not provided, and \code{x} is of length greater than one, the range of values of \code{x} is used to
#'  decide the units displayed. Numerics are also accepted - they are coerced to GRtime.
#'
#' \code{axis.GRtime} is used to plot GRtime axis, choosing, by default, breakpoints that give 'pretty' sub
#' intervals. Note that \code{\link{plot.default}} uses \code{axis.GRtime} by default if supplied with a
#' GRtime object in one of the directions. However, \code{\link{image.default}} based functions do not use
#' the class axis functions, so axes must be plotted manually.
#'
#' \code{pretty.GRtime} computes 'pretty' breakpoints, using the algorithm of \code{pretty.POSIXt}.
#' Attributes are preserved.
#'
#' \itemize{
#' \item For \code{convert.time}, \code{as.GRtime} and \code{pretty.GRtime}, a GRtime object.
#' \item For \code{format.GRtime} a character string representation.
#' \item For \code{axis.GRtime} a list containing positions and labels for axis markers.
#' }
#'
#' @importFrom graphics par axis
#'
#' @seealso \code{\link{parse.time}}, \code{\link{get.intervals}}, \code{\link{AccData}}
#'
#' @examples
#' as.GRtime("00:01")
#' #format is automatically set
#' convert.time(1:10)
#' convert.time(1:10*1000)
#' #we add a different default format
#' convert.time(1:10*1000, "\%H:\%M:\%OS3") -> t
#' t
#' str(t)
#' # we override format with our own
#' format(t, format = "\%a \%d/\%m/\%y \%H:\%M:\%OS3")
#'
#' # plot calls axis.GRtime automatically. Notice
#' # that the format attribute is used.
#' plot(t, 1:10)
#' #strip out the default format
#' t2 = convert.time(t, format = NULL)
#' plot(t2, 1:10)
#'
#' # image plots are a bit more complex
#'
#' Z = matrix(rnorm(100), 10)
#' image(x = t, y = t2, z = Z, axes = FALSE)
#' Axis(x = t, side = 1) #Axis also works
#' box() #complete the bounding box
#'
#' # custom axes
#' plot(t2, 1:10, xaxt = "n")
#'
#' @export
#'

as.GRtime <- function(x, format = NULL, ...){
  if (is.character(x)){
    return(convert.time(parse.time(x, ...), format))
  }
  else {
    # try to coerce
    return( convert.time(as.numeric(x, ...), format))}
}

#' @export

c.GRtime <- function(..., recursive = FALSE){
  structure(c(unlist(lapply(list(...), unclass))), class = "GRtime")
}

#' @export

Summary.GRtime <- function(x,..., na.rm){
  ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
  att = attributes(x)
  if (!ok)
    stop(.Generic, " not defined for GRtime objects")
  val <- NextMethod(.Generic)
  attr(val, "format") = att$format
  class(val) <- att$class #]oldClass(list(...)[[1L]])
  val
}

#' @export

Ops.GRtime <-function (e1, e2){
  if (nargs() == 1)
    stop("unary ", .Generic, " not defined for GRtime objects")
  boolean <- switch(.Generic, `<` = , `>` = , `==` = , `!=` = ,
                    `<=` = , `>=` = TRUE, FALSE)
  if (boolean) {
    if (is.character(e1))
      e1 <- parse.time(e1)
    if (is.character(e2))
      e2 <- parse.time(e2)
    return(NextMethod(.Generic))
  }

  arith <- switch(.Generic, `+` = , `-` = , `*` = , `/` = , `%/%`=
                    ,  `%%` = TRUE, FALSE)
  unitpres <- switch(.Generic, `+` = , `-` = ,  `%%` = TRUE, FALSE)

  if (inherits(e1, "GRtime")){
    att = attributes(e1)
    if ((!unitpres) && inherits(e2, "GRtime")){ att = NULL}
  } else {
    att = attributes(e2)
  }
  if (!arith) stop(.Generic, " not defined for GRtime objects")

  val = NextMethod(.Generic)
  attributes(val) = att
  return(val)
}

#' @export

`[.GRtime` <- function (x, ..., drop = TRUE)
{
  cl <- oldClass(x)
  class(x) <- NULL
  val <- NextMethod("[")
  class(val) <- cl
  val
}

#' @export

print.GRtime <- function(x, quote = FALSE, format, ...){
  if(!as.logical(length(x))) {
    cat("GRtimes(0)\n")
    return(invisible(x))
  }
  if (missing(format)){ format = NULL}
  xo <- x
  x <- format.GRtime(x, format = format)
  print(x, quote = quote)
  invisible(xo)
}

#' @export

format.GRtime <- function(x,
                            format = NULL,...){
  #format. = "h:m:s", simplify = FALSE, ...)
  if(!as.logical(length(x)))
    return("")
  if(all(is.na(x)))
    return(rep("NA", length = length(x)))
  if(!is.numeric(x))
    stop(paste(deparse(substitute(x)), "must be numeric"))
  att <- attributes(x)
  if (is.null(format)) format = att$format #maybe x has a preset format?
  if (is.null(format)){
    #choose a smart format depending on the time range?
    #most things, or single times: h:m
    #short interval < 10 minutes h:m:s
    #supershort < 10 seconds m:s:ms
    #long >24 hours d h:m
    #very long >  72 hours
    format = "%H:%M"
    if (length(x) > 1){
      rng = diff(range(x, na.rm = TRUE))
      if (rng < 10){
        format = "%M:%OS3"
      } else if (rng < 10*60) {
        format = "%H:%M:%S"
      } else if (rng < 24*60*60){
        format = "%H:%M"
      } else if (rng < 7*60*60*24){
        format = "%a %H:%M"
      } else {
        format = "%d/%m %H:%M"
      }
    }
  }

  att$class <- att$format <- NULL
  ## <NOTE>
  ## DJ's design is that
  ##   times greater than 1 day  should format like numerics
  ## To change this (e.g., have times(1.5) format as 36:00:00), simply
  ## comment the code below, and make the corresponding change in
  ## print.times().
  out = format(as.POSIXct(as.numeric(x), origin = "1970-1-1", tz = "UTC"), format, ...)

  out[x == Inf] <- "Inf"
  out[x ==  - Inf] <- "-Inf"
  attributes(out) <- att
  out
}

#' @export

pretty.GRtime <- function(x, n = 5, ...) {
  att = attributes(x)
  attributes(x) = NULL
  x = as.numeric(pretty(as.POSIXct(x, origin = "1970-1-1", tz = "UTC"), n,...))
  attributes(x) = att
  return(x)
}

#' @method axis GRtime
#'
#' @export

axis.GRtime <-function(side,
                       x = NULL,
                       at = NULL,
                       format = NULL,
                       labels  = TRUE,
                       add = TRUE,  ...){
  if (is.null(at)){
    att <- NULL
    if (inherits(x, "GRtime")) att <- attributes(x)
    bad <- (is.na(x) | abs(as.vector(x)) == Inf)
    if(side == 1 || side == 3){
      rng =  par("usr")[1:2]
      n = par("xaxp")[3] +1}
    else {
      rng = par("usr")[3:4]
      n = par("yaxp")[3] +1}
    tmp <- c(rng, as.numeric(x[!bad]))
    attributes(tmp) = att
    at = pretty.GRtime(tmp, n)}
  else {
    if(!inherits(at, "GRtime")) at <- convert.time(at)
  }
  if(missing(labels) || (is.logical(labels) && labels)){
    labels <- format(at, format = format)}
  if(add){
    axis(side, at = at, labels = labels, ...)}
  invisible(list(side = side, at = at, labels = labels))
}
