#' @name get.intervals
#'
#' @aliases get.intervals
#' @aliases print.VirtAccData
#' @aliases VirtAccData
#'
#' @title Extract an interval of data.
#'
#' @description Function for extracting sub intervals of data, and implementation of just-in-time loading.
#'
#' @usage get.intervals(x, start=0, end = 1, length = NULL,
#' time.format = c("auto", "seconds", "days", "proportion", "measurements", "time"),
#' incl.date = FALSE, simplify = TRUE ,read.from.file=FALSE, size=Inf, ...)
#'
#' @param x Object to process. Can be array,
#' @param start Start of interval.
#' @param end End of interval.
#' @param length Length of interval.
#' @param time.format Method with which \code{start} and \code{end} should be understood.
#' @param incl.date logical. Include a column denoting time?
#' @param simplify logical. If TRUE, output an array. Otherwise output a AccData object.
#' @param read.from.file logical. If TRUE, re-read the relevant time interval from the original bin file.
#' @param size {Desired number of samples in output.}
#' @param ... Additional arguments to be passed to \code{\link{read.bin}}, if \code{read.from.file} is TRUE.
#'
#' @details The function extracts the desired analysis time window specified by \code{start} and \code{end}.
#' If length is specified, then the end is set to a point \code{length} units after start.
#' The times are interpreted in terms of \code{time.format}. For convenience, a variety of time
#' window formats are accepted: \itemize{
#' \item "seconds": Seconds since start of dataset.
#' \item "days": Days since start of dataset.
#' \item "proportion": Proportional point within dataset, given as a numeric between 0 and 1.
#' \item "measurements": Raw number of samples since start of dataset.
#' \item "time": Time string, as understood via \code{\link{parse.time}}.
#' \item "auto": Default - attempt to determine time format from size and type of \code{start}.
#' }
#'
#' Some capacity for using mixed types of inputs for \code{start} and \code{length} in particular is present.
#'
#' The input object \code{x} is typically an "AccData" object, though arrays are also accepted. "VirtAccData"
#' are dealt with by using the timestamp and call information recorded within them to do a new read of the
#' original bin file, assuming this is still available. This is useful for 'just in time' reads of data.
#' "AccData" can be dealt with in this way by setting \code{read.from.file}.
#'
#' Note that for \code{read.from.file}, only "time" and "proportion" \code{time.format} are presently supported.
#'
#' With \code{simplify = FALSE}, an "AccData" S3 object with the desired records.
#'  Otherwise, an array containing either 3 or 4 columns, containing the x, y, z acceleration vectors and optionally a time vector.
#'
#' @seealso \code{\link{read.bin}}, \code{\link{AccData}}, \code{\link{get.intervals}}
#' @examples
#'
#' binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]
#'
#' #Read in a highly downsampled version of the file
#' procfile<-read.bin(binfile, downsample = 100)
#' print(procfile)
#'
#' #Overlay some segments in different colour
#' lines(get.intervals(procfile, start = 0.4, end = 0.5,
#'                     time.format = "prop", incl.date = TRUE)[,1:2],
#'                     col=2)
#'
#' lines(get.intervals(procfile, start = 0.4, end = 5,
#'                     time.format = "sec", incl.date = TRUE)[,1:2],
#'                     col=3)
#'
#' lines(get.intervals(procfile, start = "16:51", end = "16:52",
#'                     time.format = "time", incl.date = TRUE)[,1:2],
#'                     col=4)
#'
#' # Note that measurements will depend on the downsampling rate,
#' # not the original sampling rate of the data
#' lines(get.intervals(procfile, start = 100, length = 10,
#'                     time.format = "measurement", incl.date = TRUE)[,1:2],
#'                     col=5)
#'
#' #This is also understood
#' lines(get.intervals(procfile, start = "16:52:10", 30,
#'                     incl.date = TRUE)[,1:2],
#'                     col=6)
#'
#' #Now load in virtually
#' virtfile<-read.bin(binfile, virtual = TRUE)
#' #Notice that get.intervals with simplify = FALSE gives a genuine AccData object
#' realfile = get.intervals(virtfile, start = 0.5, end = 1, simplify = FALSE)
#' virtfile
#' realfile
#' #get.intervals calls read.bin automatically
#' points(get.intervals(virtfile, start = "16:52:10", "16:52:40",
#'                      incl.date = TRUE)[,1:2], col=4, pch = ".")
#'
#' #Alternatively, re-read procfile at a different resampling rate.
#' lines(get.intervals(procfile, start = "16:49:00", "16:49:30",
#'                     incl.date = TRUE, read.from.file = TRUE, downsample = 300)[,1:2],
#'                     col=2)
#'
#' @export

get.intervals = function(x,
                         start = 0,
                         end = 1,
                         length = NULL,
                         time.format = c("auto", "seconds", "days", "proportion", "measurements", "time"),
                         incl.date = FALSE,
                         simplify = TRUE,
                         read.from.file = FALSE,
                         size = Inf, ...){

  if (inherits(x, "VirtAccData")){
    read.from.file = TRUE
  }
  #virtual database, go get the relevant period first

  if (read.from.file){
    #rerun read.bin to grab the relevant part of the data
    #note: only time-like or proportion addresses of start and end work at this point....
    argorig = x$call
    readargs = c("gain", "offset", "luxv", "voltv","warn", "verbose", "do.temp", "calibrate", "downsample", "blocksize")

    argl<-as.list(match.call())
    argind<-pmatch(names(argl),readargs)
    readargs = readargs[na.omit(argind)]
    argind<-which(!is.na(argind))
    argorig = c( argl[argind], argorig)
    argorig = argorig[which((!duplicated(names(argorig))) & (names(argorig)!= ""))]
    argorig$virtual = F
    argorig$start= start
    argorig$end = end

    x = do.call(read.bin, args = argorig)

    start = 0
    end = 1
    time.format = "proportion"
  }

  sampling.freq = 100
  time.format = match.arg(time.format)


  if (length(start) > 1) {
    end = (start)[2]
    start = (start)[1]
  }
  #auto detect time format
  if (time.format == "auto"){
    if (is.character(start)){
      time.format = "time"
    }
    else if (start <1){
      time.format = "proportion"
    }
    else if (floor(start) == start) {
      time.format = "seconds"
    }
    else {
      time.format = "days"
    }
  }

  if (is.list(x)){
    if ((time.format == "time")||(time.format =="seconds")){
      times = x[,1]
    }

    sampling.freq = x$freq
    if (simplify) {
      x = x$data.out[,(2- incl.date):4]
    }
  } else {
    if (ncol(x) == 3){
      x = cbind (1:nrow(x)/sampling.freq, x)
    }
    if ((time.format == "time") ){
      times = x[,1]
    }
    if (!incl.date){
      x = x[,-1]
    }
  }

  if (simplify){
    n = nrow(x)
  } else {
    n = nrow(x$data.out)
  }

  if (time.format == "time"){
    #if (nchar(start) >= 11){
    #tmp = strsplit(start, " ")
    #tmp = tmp[[1]][which( nchar(tmp[[1]]) >= 5)]
    #if (length(tmp) > 1){
    #start = tmp[1]
    #end = tmp[2]
    #}
    #}
    #todo: make this work

    start = parse.time(start, format = "seconds")
    t1midnight = floor(times[1] / (60*60*24)) * 60*60*24
    t1 = times[1]

    if (start < t1midnight){
      start = start + t1midnight
    }

    if (start < t1){
      start = start + 60*60*24
    }

    start = findInterval ( start, times, all.inside = T)
    t1 = times[start+1]

    if (is.character(end)){
      end = parse.time(end, format = "seconds")
      if (end < t1midnight){
        if (end >= 24*60*60){
          end = end + t1midnight
        }
        else {
          end = end +ceiling((t1 - end)/(60*60*24)) * 60*60*24
        }
      }
      end = findInterval(end, times, all.inside = T) + 1
    }
    else {
      length = end * sampling.freq
      end = NULL
    }

    time.format = "measurements"

  }

  #if (time.format == "date"){
  #start = (start - times[1]) * 60*60*24
  #if (inherits(end, "times2")){
  #end =(end - times[1]) * 60*60*24
  #}
  #time.format = "seconds"
  #}


  if (is.null(length)){
    if ( end < start){
      length = end
      end = NULL
    }
  }

  if (!is.null(length)) {
    if ((time.format == "proportion") && (length >= 1)){
      time.format ="seconds"
      start = (start * n/sampling.freq)
    }
    end = start + length
  }

  #convert into measurements

  if (time.format == "proportion"){
    start = ceiling(start * n)
    end  = floor(end * n)
  } else if (time.format == "seconds") {
    if (exists("times") && (start > times[1])){
      start = findInterval(start, times)
      end = findInterval(end+0.01, times)
    } else {
      start = ceiling(start * sampling.freq)
      end = floor(end * sampling.freq)
    }
  } else if (time.format == "days"){
    start = ceiling(start * sampling.freq*60*60*24)
    end = floor(end * sampling.freq*60*60*24)
  }

  start = max(start,1)
  end = min(end, n)

  if (incl.date) cat("Extracting time interval: ", format(convert.time(c(x[start,1], x[end,1]))) , "\n")
  ind = start:end
  tmp = 1
  if (length(ind) > size){
    tmp = ceiling(length(ind) / size)
    ind = ind[(ind %% tmp == 0)]
  }

  x = x[ind,]
  if (is.list(x)) x$freq = x$freq /tmp

  return(x)
}
