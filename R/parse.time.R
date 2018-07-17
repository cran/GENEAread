#' @name parse.time
#'
#' @aliases parse.time
#'
#' @title Parses a character time representation to another format.
#'
#' @description Converts a character vector in a variety of forms into either the raw second, second classed as POSIXct, or days since Unix epoch.
#'
#'
#' @usage  parse.time(t="",format=c("seconds", "days", "POSIX"), tzone = 0,
#'  start = NULL, startmidnight = NULL)
#'
#' @param t A character string representation of a date-time expression.
#' @param format A character string indicating which representation to output.  Can be either \code{seconds}, \code{days} or \code{POSIX}.
#' @param tzone The time zone the time is given in, expressed as an offset from UTC in hours.
#' @param start Earliest allowable time stamp in the data, as seconds since Unix epoch.
#' @param startmidnight Midnight of day '0' in the data, as seconds since Unix epoch.
#'
#' @details
#' The function processes character vectors of the form "DATE TIME" -- that is to say, a maximum of two terms separated by a space per value.
#'
#' "TIME" is given in 24 hour format, seperated by colons, in "hh:mm", "hh:mm:ss", "hh:mm:ss:ms" or "hh:mm:ss.ms" format. If ommitted, the time is taken to be 00:00:00.000.
#'
#'  "DATE" can be a date representation as "YYYY-MM-DD", "DD/MM/YY" or "DD/MM/YYYY" (noting the use of a colon or backslash seperator to distinguish between the two). Alternatively, with \code{start} and/or \code{startmidnight} supplied, an integer "NN" or string "DOW" corresponding to a day of the week can be used instead. Then, the function will find the first timestamp matching the correct "TIME", that falls NN midnights after \code{startmidnight} and is after \code{start}, or, in the latter case, the first timestamp after the day of \code{start} that matches the appropriate day of the week. If a blank "DATE" is supplied, the function will either use the UNIX epoch, or find the first match, corresponding to the case NN = 0.
#'
#'  Once this is done the time is converted to the required format: \code{POSIX} is the usual R POSIXct format; \code{days} is the julian days since UNIX epoch 1970-1-1; \code{seconds} is the number of seconds (including subseconds) since 1970-1-1. Note that for formats other than POSIX, the output is in the same timezone as \code{tzone}. POSIX stores the time internally as the time in UTC, and applies a format that gives this time local to the user.
#'
#' @return A converted date-time string in the specified format. In the case of "seconds", or "days", a numeric. For POSIX, a \code{\link{POSIXct}} object.
#'
#' @seealso \code{\link{convert.time}}, \code{\link{get.intervals}}
#'
#' @examples
#' t1 = parse.time("2012-06-21 13:04:01"); print(t1)
#' parse.time("21/06/12 13:04:01") #gives the same result
#'
#' parse.time(c("19/07/70", "20/07/70"), format = "days")
#' #results here will depend on your locale
#' parse.time(c("19/07/70", "20/07/70"), format = "POSIX", tzone = -4)
#'
#' #one is the same day, one can only find a match the next day
#'  parse.time("13:05", start = t1) - t1
#'  parse.time("13:00", start = t1) - t1
#'  #asking to wait 1 midnight means both times are considered as
#'  #times on the same, full day of data
#'  parse.time(c("1 13:05", "1 13:00"), start = t1) - t1
#'  #2012-06-21 is a Thursday, so this is equivalent
#'  parse.time(c("Fri 13:05", "Fri 13:00"), start = t1) - t1
#'  #Longer form days of the week are also understood. Note that
#'  #the first day does not get matched.
#'  parse.time(c("Thursday 13:05", "Thursday 13:00"), start = t1) - t1
#'
#' @export

parse.time <- function(t="",
                       format = c("seconds", "days", "POSIX"),
                       tzone = 0,
                       start = NULL,
                       startmidnight = NULL){

  dow = NULL
  format = match.arg(format)
  millisec = rep(0, length(t))
  offset = 0
  t1= t[1]
  informat = ""

  #do we have time in here?
  switch(
    length(strsplit(t1, split = ":")[[1]]),

    "1" = {
      informat = ""
      },

    "2" = {
      informat = "%H:%M"
      },

    "3" = {
      informat = "%H:%M:%S"
      if (length(strsplit(t1, split = "\\.")[[1]])>1){
        millisec = sapply( strsplit(t, split="\\."), function(t) as.numeric( paste("0.", t[2], sep = "")))
        t = sapply( strsplit(t, split="\\."), function(t) t[1])
        }
      },

    "4" = {
      informat = "%H:%M:%S"
      millisec = sapply( strsplit(t, split=":"), function(t) as.numeric( paste("0.", t[4], sep = "")))
      t = sapply( strsplit(t, split=":"), function(t) paste(t[1:3], collapse = ":"))
      }
    )

  #do we have date?
  #strip whitespace/time
  t1 = strsplit(t1, split=" ")[[1]]
  t1 = t1[min(which(t1 != ""))]
  #mode yyyy-mm-dd
  if (length(strsplit(t1, split = "-")[[1]]) > 1){
  	informat = paste("%Y-%m-%d", informat, sep = " ")
  } else if (length(strsplit(t1, split = "/")[[1]]) > 1){
  #mode d/m/y
  	if (nchar(strsplit(t1, split = "/")[[1]][3])<4){
  		informat = paste("%d/%m/%y", informat, sep = " ")
  	} else {
  		informat = paste("%d/%m/%Y", informat, sep = " ")
  	}
  } else if (length(strsplit(t1, split = ":")[[1]]) > 1) {
  #nothing found, first non-open field is time
  #send us back to 1970
  	offset = as.numeric(strptime("00:00", format = "%H:%M", tz = "UTC"))
  } else {
  #add some days?
  	t = lapply(strsplit(t, split=" "), function(t) (t[t != ""]))
  	if (!suppressWarnings(is.na(as.numeric(t[[1]][1])))){
  		t1 = sapply(t, function(x) as.numeric(x[ 1]))
  		t = sapply(t, function(t) t[ 2])
  	} else {
  		dow = sapply(t, function(x) x[1])
  		t = sapply(t, function(t) t[ 2])
  		t1 = 0
  	}
  		offset = as.numeric(strptime("00:00", format = "%H:%M", tz = "UTC")) - t1 * 24*60*60
  }

  if (informat != ""){
  	t = as.POSIXct(strptime(t, format = informat), tz = "UTC")
  } else {
  	t = as.POSIXct(strptime("00:00", format = "%H:%M", tz = "UTC"))
  }

  t= t+ millisec - offset

  if ((!is.null(start)) || (!is.null(startmidnight))){
  	if (is.null(startmidnight)) startmidnight = floor(start/(60*60*24)) * 60*60*24
  	if (is.null(start)) start = startmidnight
  #resolve ambiguity
  	if (is.null(dow)){
  		if (t[1] < startmidnight) {
  		  t= t + startmidnight
  		}
  		if (t[1] < start){
  		  t= t + ceiling((start-as.numeric(t[1]))/(60*60*24)) * 60*60*24
  		}
  	} else {
  		#day of the week processing
  		#get DOW for midnight on start
  		startmidnight= as.POSIXct(floor(start/(60*60*24)) * 60*60*24,
  		                          origin = "1970-1-1",
  		                          tz = "UTC")
  		next_week <- as.Date(startmidnight) + 1:7
  		dow = substr(tolower(dow),
  		             1,
  		             min(nchar(dow)))
  		startmidnight = as.numeric(startmidnight)
  		t = t + startmidnight + (match(dow,
  		                               substr(tolower(weekdays(next_week)),
  		                                      1,
  		                                      min(nchar(dow)))) - 0) * 24*60*60

  	}
  }

  if (format == "seconds"){
  	t = as.numeric(t)
  } else if (format == "POSIX"){
  	t = c(t - tzone*60*60)
  } else {
  	t = as.numeric(t) / (60*60*24)
  }
  return(t)
}


