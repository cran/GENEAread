
#' @name AccData
#' @title Accelerometer Data Onject
#' @description Accelerometer Data Output from read.bin function
#' @docType data
#' @format An AccData object
#' @source Output of \code{\link{read.bin}}
#' @keywords datasets
#' @seealso \code{read.bin}
#' @examples
#' binfile  = system.file("binfile/TESTfile.bin", package = "GENEAread")[1]
#' #Read in the entire file, calibrated
#' procfile<-read.bin(binfile)
#' print(procfile)
#'
#' plot(procfile$temperature,
#' xlim = c(min(procfile$data.out[,1]),
#'         max(procfile$data.out[,1])),
#' ylim = c(10,40))
#'
#' plot(procfile$data.out[,c(1,7)])

setClass("AccData", representation(data.out = "matrix",
                                   page.timestamps = "numeric",
                                   freq = "numeric",
                                   filename = "character",
                                   page.numbers = "numeric",
                                   call = "list",
                                   page.volts = "numeric",
                                   pagerefs = "numeric",
                                   header = "list"))

print.AccData <- function(x){
  cat("GENEAread dataset: ", nrow(x$data.out), "records at", round(x$freq,2),
      "Hz (Approx ", round(object.size(x$data.out)/1000000) ,"MB of RAM)\n")

  cat(format.GRtime(x$data.out[1,1], format = "%y-%m-%d %H:%M:%OS3 (%a)")," to ",
      format.GRtime(tail(x$data.out[,1],1),format = "%y-%m-%d %H:%M:%OS3 (%a)"), "\n")
  #}
  cat("[", x$filename, "]\n")
}

setMethod("print.AccData", signature(x = "AccData"), function(x){
  cat("GENEAread dataset: ", nrow(x$data.out), "records at", round(x$freq,2),
      "Hz (Approx ", round(object.size(x$data.out)/1000000) ,"MB of RAM)\n")

  cat(format.GRtime(x$data.out[1,1], format = "%y-%m-%d %H:%M:%OS3 (%a)")," to ",
      format.GRtime(tail(x$data.out[,1],1),format = "%y-%m-%d %H:%M:%OS3 (%a)"), "\n")
  #}
  cat("[", x$filename, "]\n")
})

"[.AccData"    <- function (x,
                            i = 1:dim(x$data.out)[1],
                            j = NULL,
                            drop = TRUE) {
  if (is.null(j)){
    x$page.timestamps = x$page.timestamps[ unique(ceiling(i/300))]
    x$data.out = x$data.out[i,]
    return(x)
  }
  if ((length(j) == ncol(x$data.out) )&& (max(j) <= 1)){
    j = which(j)
  }
  if ( j[1] == 1 ){

    if (length(j) != 1){
      value = x$data.out[i, j[-1] , drop = F]

      return( data.frame( time = convert.time(x$data.out[i,1, drop = T]), value  ))
    } else{
      return (convert.time(x$data.out[i,1, drop = drop]))
    }
  } else {
    return(x$data.out[i,j, drop=drop])
  }
}

setMethod("[.AccData", signature(x = "AccData"), function (x,
                                                           i = 1:dim(x$data.out)[1],
                                                           j = NULL,
                                                           drop = TRUE) {
  if (is.null(j)){
    x$page.timestamps = x$page.timestamps[ unique(ceiling(i/300))]
    x$data.out = x$data.out[i,]
    return(x)
  }
  if ((length(j) == ncol(x$data.out) )&& (max(j) <= 1)){
    j = which(j)
  }
  if ( j[1] == 1 ){

    if (length(j) != 1){
      value = x$data.out[i, j[-1] , drop = F]

      return( data.frame( time = convert.time(x$data.out[i,1, drop = T]), value  ))
    } else{
      return (convert.time(x$data.out[i,1, drop = drop]))
    }
  } else {
    return(x$data.out[i,j, drop=drop])
  }
})

"$.AccData" <- function(x, name){
  nmatch <- try(match.arg(name,
                          c("time", "x", "y", "z", "xyz",
                            "temperature", "button", "voltage",
                            "light", "svm")),
                silent = TRUE)
  if (inherits(nmatch, "try-error")){
    class(x) <- NULL
    return(x[[name, exact = FALSE]])
  } else {
    #	x = unclass(x)
    ind = switch(nmatch, time = 1, x = 2, y = 3, z = 4, xyz = 2:4,
                 temperature = 7, button = 6, light = 5, voltage = 8, svm = 9)
    if (identical(ind, 8)){
      return(rep(x$page.volt,
                 each = ceiling(nrow(x)/length(x$page.volt)) )[1:nrow(x)])
    } else if (identical(ind, 9)){
      return(svm(x))
    } else {
      return(x[,ind])
    }
  }
}

setMethod("$.AccData", signature(x = "AccData"), function(x, name){
  nmatch <- try(match.arg(name,
                          c("time", "x", "y", "z", "xyz",
                            "temperature", "button", "voltage",
                            "light", "svm")),
                silent = TRUE)
  if (inherits(nmatch, "try-error")){
    class(x) <- NULL
    return(x[[name, exact = FALSE]])
  } else {
    #	x = unclass(x)
    ind = switch(nmatch, time = 1, x = 2, y = 3, z = 4, xyz = 2:4,
                 temperature = 7, button = 6, light = 5, voltage = 8, svm = 9)
    if (identical(ind, 8)){
      return(rep(x$page.volt,
                 each = ceiling(nrow(x)/length(x$page.volt)) )[1:nrow(x)])
    } else if (identical(ind, 9)){
      return(svm(x))
    } else {
      return(x[,ind])
    }
  }
})

c.AccData <- function(x, ...){
  tmp = list(x)
  out = list()
  out$data.out = NULL
  out$page.timestamps = NULL
  for (i in 1:length(tmp)){
    out$data.out = rbind(out$data.out, tmp[[i]]$data.out)
    out$page.timestamps = c(out$page.timestamps, tmp[[i]]$page.timestamps)
  }
  out$freq = tmp[[1]]$freq
  out$filename = tmp[[1]]$filename
  class(out) = class(tmp[[1]])
  return(out)
}

setMethod("c.AccData", signature(x = "AccData"), function(x, ...){
  tmp = list(...)
  out = list()
  out$data.out = NULL
  out$page.timestamps = NULL
  for (i in 1:length(tmp)){
    out$data.out = rbind(out$data.out, tmp[[i]]$data.out)
    out$page.timestamps = c(out$page.timestamps, tmp[[i]]$page.timestamps)
  }
  out$freq = tmp[[1]]$freq
  out$filename = tmp[[1]]$filename
  class(out) = class(tmp[[1]])
  return(out)
})
