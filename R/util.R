
#' @title  Constrain
#'
#' @description constrains a dataset
#'
#' @param x dataset that is to be constrained
#' @param minimum Value of the lowest acceptable number, before being replaced by this lower bound
#' @param maximum Value of the highest acceptable number, before being replaced by this upper bound
#'
#' @keywords internal

constrain <- function(x, minimum, maximum){
  if (missing(maximum)) maximum = Inf
  if (missing(minimum)) minimum = -Inf

  if (minimum > maximum) {
    temp <- minimum
    minimum <- maximum
    maximum <- temp
  }
  x = replace(x, which(x > maximum), maximum)
  x = replace(x, which(x < minimum), minimum)
  x
}

#' @title shift
#'
#' @description shift an object x in a direction/expansion
#'
#' @param x data object to process
#' @param offset 2 integers to offset a dataset
#' @param expand 2 integer list to expand the dataset
#' @param fill method to run on the shift
#'
#' @keywords internal

shift = function(x, offset, expand = 0, fill = c("zero", "edge", "loop", "mean")){
  fill = match.arg(fill)

  if (length(dim(x)) == 2){
    if (length(expand) == 1) expand = rep(expand,2)
    n = nrow(x)
    p = ncol(x)
    xlim = nrow(x) + expand[1]
    ylim = ncol(x) + expand[2]
    x2 = matrix(0, xlim, ylim)
    if (fill == "mean"){
      x2 = matrix(mean(x), xlim, ylim)
    }
    TL = c(1,1) - pmin(0, offset)
    #BR = pmin(c(nrow(x) , ncol(x)) - TL + offset + c(1,1), c(xlim, ylim)) - offset
    BR = pmin(offset + c(n,p) , c(xlim, ylim)) - offset

    x2[0:((BR-TL)[1]) +max(offset[1],0)+1,0:((BR-TL)[2])+max(offset[2],0)+1] =  x[TL[1]:BR[1], TL[2]:BR[2]]

    if ((fill != "zero") && (fill != "mean") ){

      #top left corner
      #if (min(offset) > 0){
      #
      #if (fill == "loop"){
      #	x2[1: offset[1], 1:offset[2]] = x[ 1:offset[1] + n - offset[1], 1:offset[2] + p - offset[2]]
      #} else {
      #x2[1: offset[1], 1:offset[2]] = x[1,1]
      #}
      #}

      #top edge
      if (offset[1] > 0){
        x2[1:offset[1],] = matrix(shift(x[1,],
                                        offset[2],
                                        expand = ylim - p,
                                        fill="edge"),
                                  nrow  =  offset[1],
                                  byrow = T, ncol = ylim)
      }
      #left edge

      if (offset[2] > 0){
        x2[,1:offset[2]] = matrix(shift(x[,1],
                                        offset[1],
                                        expand = xlim - n,
                                        fill="edge"),
                                  ncol = 1)
      }

      #bottom edge
      if (n + offset[1] < xlim){
        x2[(n + offset[1]+1) :xlim,] = matrix(shift(x[n,],
                                                    offset[2],
                                                    expand = ylim - p,
                                                    fill="edge"),
                                              nrow = xlim - offset[1] - n,
                                              ncol= ylim, byrow = T)
      }

      #right edge

      if (p + offset[2] < ylim){
        x2[,(p + offset[2]+1) :ylim] = matrix(shift(x[,p],
                                                    offset[1],
                                                    expand = xlim - n,
                                                    fill="edge"),
                                              ncol = 1)
      }

      if (fill == "loop"){
        print("Currently umimplemented loop fill, edging instead")
        #if (sum(expand) != 0)	print("Warning, loop results may not make sense if matrix size changes")
      }
    }
  }
  else {
    #return(drop(shift(x=matrix(x, ncol=1), offset = c(offset,0), expand =c(expand, 0), fill = fill)))
    offset = offset[1]
    expand  = expand[1]
    x2 = rep(0, length(x) + expand)
    beg = 1 - pmin(offset, 0)
    end =  min (offset + length(x), length(x2)) - offset # min(  length(x) - beg + offset + 1, length(x) + expand) - offset

    x2[( max(offset,0) +1) : min(length(x2), length(x)+ offset)] = x[beg:end]

    if (fill == "loop"){
      if (offset > 0){
        x2[1:offset] = tail(x, offset)
      }
      else{
        x2[((length(x) + offset + 1): length(x2))] = x[1:(-offset)]
      }
    }
    else if (fill == "edge"){
      if (offset > 0){
        x2[1:offset] = x[1]
      }
      if ((length(x) + offset) <  length(x2)){
        x2[ (length(x) + offset +1) : length(x2)] = tail(x,1)
      }
    }
  }
  x2
}

#' @title convert vector
#'
#' @description puts a vector into the range 0-1
#'
#' @keywords internal

conv01 <- function(x){
  (x - min(x))/ (max(x)- min(x))
}

#' @title \%bq\%
#'
#' @description quantile version of bt
#'
#' @param x First object to pass
#' @param y Second object to pass
#' @importFrom stats quantile
#'
#' @keywords internal

"%bq%" = function(X, y){
  if (is.character(y)){
    if (length(y) == 4){
      y = paste(y[1],y[2], ",",y[3], y[4], sep="")
    }
    y = strsplit(y, ",")[[1]]
    nc = nchar(y[2])
    yl = quantile(X, c(as.numeric(substring(y[1], 2)), as.numeric(substring(y[2], 1,nc - 1))))
    if (substr(y[1],1,1) == "["){
      res = (X >= yl[1])
    }else {
      res = (X >  yl[1])
    }
    if (substr(y[2],nc,nc) == "]"){
      res = res &(X <= yl[2] )
    }else {
      res = res & (X < yl[2])
    }
  } else {
    y = quantile(X, y)
    res = (X >= y[1] ) & (X<= y[2])
  }
  res
}


#' @title \%bt\%
#'
#' @description between' operator for convenience
#' takes \code{[min, max)}, or \code{c("[", min, max, "]")} style second terms
#' default is \code{[min, max]} for c(,) terms
#'
#' @param x First object to pass
#' @param y Second object to pass
#'
#' @keywords internal

"%bt%" = function(X, y){
  if (is.character(y)){
    if (length(y) == 4) y = paste(y[1],y[2], ",",y[3], y[4], sep="")
      y = strsplit(y, ",")[[1]]
      if (substr(y[1],1,1) == "["){
        res = (X >= as.numeric(substring(y[1], 2)))
      } else {
      res = (X > as.numeric(substring(y[1], 2)))
      }
    nc = nchar(y[2])
    if (substr(y[2],nc,nc) == "]"){
      res = res &(X <= as.numeric(substring(y[2],1, nc -1)))
    } else {
      res = res & (X < as.numeric(substring(y[2], 1,nc - 1)))
    }
  } else {
    res = (X >= y[1] ) & (X<= y[2])
  }
  res
}

#' @title removezero
#'
#' @description removes zero's from an object
#'
#' @param obj object to remove zeros from
#'
#' @keywords internal

removeZero <- function(obj){
  obj[which(obj!=0)]
}

#' @title seq.log
#'
#' @description Create a sequence of log values
#'
#' @param from Value to Create the sequence from. default to 1.
#' @param to Value to create the sequence to. default to 1
#' @param length.out length of list outputted
#' @param add.zero include 0 at the start
#' @param shifting Shifting the sequence by this value
#'
#' @keywords internal

seq.log <- function(from = 1,
                    to = 1,
                    length.out = 50,
                    add.zero = FALSE,
                    shifting = 0){

  res = exp(seq(from = log(from + shifting),
                to = log(to + shifting),
                length=length.out - add.zero)) - shifting
  if (add.zero) {
    if (from > to) {
      res = c(res,0)
    } else {
      res = c(0,res)
    }
  }
  res
}

#' @title expand
#'
#' @description expand an obecjt x
#'
#' @param x object to expand
#' @param length length to expand object
#'
#' @importFrom utils tail
#'
#' @keywords internal

expand <- function(X, length = (length(X)*100)){
  c(rep(X, each = floor(length / length(X))), rep(tail(X,1), length - length(X) * floor(length/length(X))))
}

#' @title pseudoline
#' @description  Plot a line graph, with breaks when things change too much
#' (assume x is sorted)
#'
#' @param x sorted vector
#' @param y response variable
#' @param max.shift maximum to shift data
#' @param new Create new plot
#'
#' @keywords internal


pseudolines <- function(x, y = NULL, max.shift, new = FALSE,...){
  if (new){
    plot(x,y, type = "n", ...)
  }
  if (is.null(y)){
    if ((length(dim(x)) > 1) && (ncol(x)==2 ) ){
      y = x[,2]
      x = x[,1]
    }
    else {
      y = drop(x)
      x = 1:length(y)
    }
  }

  n = length(y)

  if (missing(max.shift)){
    max.shift = max(sqrt(2) /n, 0.05)
  }

  xrng = par("usr")
  yrng = xrng[4] - xrng[3]
  xrng = xrng[2] - xrng[1]

  deltas = sqrt(rowSums(cbind(diff(x/xrng), diff(y/yrng))^2))
  deltas = (deltas > max.shift)
  switches = c(which(diff(deltas) != 0), n)

  pos = 1
  for (i in 1: (length(switches)+1)){
    type = deltas[pos]
    leng = switches[i] - pos + 2

    if (type == TRUE){
      leng = leng - 2
      if (leng > 0){
        points(x[pos + 1:leng ], y[pos + 1:leng], ...)
      }
      leng = leng + 2
    }
    else{
      lines(x[pos -1+ 1:leng], y[pos -1+ 1:leng], ...)
    }
    pos = pos + leng-1
    if (pos > n ) {
      break
    }
  }
}

