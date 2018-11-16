
#' @title GENEActiv.calibrate
#'
#' @name GENEActiv.calibrate
#'
#' @description Function starts by identifying ten second windows of non-movement.
#' Next, the average acceleration per axis per window is used to estimate calibration error
#' (offset and scaling) per axis. The function provides recommended correction factors to address
#' the calibration error and a summary of the callibration procedure.
#'
#' @param binfile A filename of a file to process.
#' @param use.temp use temperature sensor data for calibration
#' @param spherecrit The minimum required acceleration value (in g) on both sides of 0 g for each axis.
#' Used to judge whether the sphere is sufficiently populated
#' @param minloadcrit The minimum number of hours the code needs to read for the autocalibration procedure
#' to be effective (only sensitive to multitudes of 12 hrs, other values will be ceiled).
#' After loading these hours only extra data is loaded if calibration error has not been reduced to under 0.01 g
#' @param printsummary if TRUE will print a summary when done
#' @param chunksize number between 0.2 and 1 to specificy the size of chunks to be loaded as a fraction
#' of a 12 hour period, e.g. 0.5 equals 6 hour chunks. The default is 1 (12 hrs).
#' For machines with less than 4Gb of RAM memory a value below 1 is recommended.
#' @param windowsizes Three values to indicate the lengths of the windows as in c(window1,window2,window3):
#' window1 is the short epoch length in seconds and by default 5 this is the time window over which
#' acceleration and angle metrics are calculated, window2 is the long epoch length in seconds for which
#' non-wear and signal clipping are defined, default 900. However, window3 is the window length of data used
#' for non-wear detection and by default 3600 seconds. So, when window3 is larger than window2 we use
#' overlapping windows, while if window2 equals window3 non-wear periods are assessed by non-overlapping
#' windows.
#'
#' @details
#'
#' @details The outputs from the function are as follows \itemize{
#' \item scale scaling correction values, e.g. c(1,1,1)
#' \item offset offset correction values, e.g. c(0,0,0)
#' \item tempoffset correction values related to temperature, e.g. c(0,0,0)
#' \item cal.error.start absolute difference between Euclidean norm during all non-movement windows and 1 g before autocalibration
#' \item cal.error.end absolute difference between Euclidean norm during all non-movement windows and 1 g after autocalibration
#' \item spheredata average, standard deviation, Euclidean norm and temperature (if available) for all ten second non-movement windows as used for the autocalibration procedure
#' \item npoints number of 10 second no-movement windows used to populate the sphere
#' \item nhoursused number of hours of measurement data scanned to find the ten second time windows with no movement
#' \item meantempcal mean temperature corresponding to the data as used for autocalibration.
#' }
#'
#' @author Vincent T van Hees <vincentvanhees@gmail.com>
#' Zhou Fang
#' Charles Sweetland <charles@Sweetland-solutions.co.uk>
#'
#' @references van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI, White T,
#' Wareham NJ, Brage S. Auto-calibration of accelerometer data for free-living physical activity assessment
#' using local gravity and temperature: an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
#'
#' @importFrom stats lm.wfit
#'
#' @export






GENEActiv.calibrate = function(binfile,
                               use.temp = TRUE,
                               spherecrit = 0.3,
                               minloadcrit = 72,
                               printsummary = TRUE,
                               chunksize = c(),
                               windowsizes = c(5,900,3600)) {

  if (length(chunksize) == 0) chunksize = 1
  if (chunksize > 1) chunksize = 1
  if (chunksize < 0.4) chunksize = 0.4
  filename = unlist(strsplit(as.character(binfile),"/"))
  filename = filename[length(filename)]

  # set parameters
  filequality = data.frame(filetooshort=FALSE,filecorrupt=FALSE,filedoesnotholdday = FALSE)
  ws4 = 10 # Epoch for recalibration, don't change
  ws2 = windowsizes[2] # Dummy variable
  ws =  windowsizes[3] # window size for assessing non-wear time (seconds)
  i = 1 # Counter to keep track of which binary block is being read
  count = 1 # Counter to keep track of the number of seconds that have been read
  LD = 2 # Dummy variable used to identify end of file and to make the process stop
  cal.error.start=cal.error.end=c()
  spheredata=c()
  tempoffset=c()
  npoints=c()
  scale = c(1,1,1)
  offset = c(0,0,0)
  bsc_cnt = 0
  bsc_qc = data.frame(time=c(),size=c())
  # Find the frequency of the file.
  options(warn=-1) #turn off warnings
  H = header.info(binfile)
  options(warn=0) #turn off warnings
  sf = as.numeric(H$Value[[2]])

  if (sf == 0) stop("File is corrupt")

  #creating matrixes for storing output
  S = matrix(0,0,4) #dummy variable needed to cope with head-tailing succeeding blocks of data
  NR = ceiling((90*10^6) / (sf*ws4)) + 1000 #NR = number of 'ws4' second rows (this is for 10 days at 80 Hz)
  meta = matrix(99999,NR,8) #for meta data
  blocksize = round((14512 * (sf/50)) * (chunksize*0.5))

  # Read file
  switchoffLD = 0 #dummy variable part of "end of loop mechanism"

  while (LD > 1) {
    P = c()
    if (i  == 1) {
      cat(paste("\nLoading block: ",i,sep=""))
    } else {
      cat(paste(" ",i,sep=""))
    }

    #try to read data blocks based on monitor type and data format
    options(warn = -1) #turn off warnings (code complains about unequal rowlengths

    blocknumber = i
    #### Replacing with read.bin - Need to think about how we re write this in. ####
    AccData = read.bin(binfile = binfile,
                       start = (blocksize*(blocknumber-1)),
                       end = (blocksize*blocknumber),
                       calibrate = TRUE,
                       do.temp = TRUE,
                       mmap.load = FALSE)

    options(warn = 0) #turn on warnings
    data = AccData$data.out
    if (length(data) > 0) { #would have been set to zero if file was corrupt or empty
      # add left over data from last time
      if (min(dim(S)) > 1) {
        data = rbind(S,data)
      }
      LD = nrow(data)
      #store data that could not be used for this block, but will be added to next block
      use = (floor(LD / (ws*sf))) * (ws*sf) #number of datapoint to use
      if (length(use) > 0) {
        if (use > 0) {
          if (use != LD) {
            S = as.matrix(data[(use+1):LD,]) #store left over
          }
          data = as.matrix(data[1:use,])
          LD = nrow(data) #redefine LD because there is less data
          ##==================================================
          dur = nrow(data)	#duration of experiment in data points
          durexp = nrow(data) / (sf*ws)	#duration of experiment in hrs
          # Initialization of variables
          Gx = as.numeric(data[,2]); Gy = as.numeric(data[,3]); Gz = as.numeric(data[,4]); temperaturre = as.numeric(data[,7])
          temperature = as.numeric(data[,7])

          temperaturecolumn = 7;
          lightcolumn = 5

          if (use.temp == TRUE) { #also ignore temperature for GENEActive if temperature values are unrealisticly high or NA
            if (length(which(is.na(mean(as.numeric(data[1:10,7]))) == T)) > 0) {
              cat("\ntemperature is NA\n")
              use.temp = FALSE
            } else if (length(which(mean(as.numeric(data[1:10,7])) > 40)) > 0) {
              cat("\ntemperature is too high\n")
              use.temp = FALSE
            }
          }

          #=============================================
          # non-integer sample frequency is a pain for deriving epoch based sd
          # however, with an epoch of 10 seconds it is an integer number of samples per epoch
          EN = sqrt(Gx^2 + Gy^2 + Gz^2)
          D1 = g.downsample(EN,sf,ws4,ws2)
          EN2 = D1$var2
          #mean acceleration
          D1 = g.downsample(Gx,sf,ws4,ws2); 	GxM2 = D1$var2
          D1 = g.downsample(Gy,sf,ws4,ws2); 	GyM2 = D1$var2
          D1 = g.downsample(Gz,sf,ws4,ws2); 	GzM2 = D1$var2
          if (use.temp == TRUE) {
            D1 = g.downsample(temperature,sf,ws4,ws2);
            TemperatureM2 = D1$var2
          }
          #sd acceleration
          dim(Gx) = c(sf*ws4,ceiling(length(Gx)/(sf*ws4))); 	GxSD2 = apply(Gx,2,sd)
          dim(Gy) = c(sf*ws4,ceiling(length(Gy)/(sf*ws4))); 	GySD2 = apply(Gy,2,sd)
          dim(Gz) = c(sf*ws4,ceiling(length(Gz)/(sf*ws4))); 	GzSD2 = apply(Gz,2,sd)
          #-----------------------------------------------------
          #expand 'out' if it is expected to be too short
          if (count > (nrow(meta) - (2.5*(3600/ws4) *24))) {
            extension = matrix(" ",((3600/ws4) *24),ncol(meta))
            meta = rbind(meta,extension)
            cat("\nvariabel meta extended\n")
          }
          #storing in output matrix
          meta[count:(count-1+length(EN2)),1] = EN2
          meta[count:(count-1+length(EN2)),2] = GxM2
          meta[count:(count-1+length(EN2)),3] = GyM2
          meta[count:(count-1+length(EN2)),4] = GzM2
          meta[count:(count-1+length(EN2)),5] = GxSD2
          meta[count:(count-1+length(EN2)),6] = GySD2
          meta[count:(count-1+length(EN2)),7] = GzSD2
          if (use.temp == TRUE) {
            meta[count:(count-1+length(EN2)),8] = TemperatureM2
          }
          count = count + length(EN2) #increasing "count": the indicator of how many seconds have been read
          rm(Gx); rm(Gy); rm(Gz)
          # reduce blocksize if memory is getting higher
          gco = gc()
          memuse = gco[2,2] #memuse in mb
          bsc_qc = rbind(bsc_qc,c(memuse,Sys.time()))
          if (memuse > 4000) {
            if (bsc_cnt < 5) {
              if ((chunksize * (0.8 ^ bsc_cnt)) > 0.2) {
                blocksize = round(blocksize * 0.8)

                bsc_cnt = bsc_cnt + 1
              }
            }
          }
        }
      }
    } else {
      LD = 0 #once LD < 1 the analysis stops, so this is a trick to stop it
      cat("\nstop reading because there is not enough data in this block\n")
    }
    spherepopulated = 0
    if (switchoffLD == 1) {
      LD = 0
    }
    meta_temp = data.frame(V = meta)
    cut = which(meta_temp[,1] == 99999)
    if (length(cut) > 0) {
      meta_temp = as.matrix(meta_temp[-cut,])
    }
    nhoursused = (nrow(meta_temp) * 10)/3600
    if (nrow(meta_temp) > (minloadcrit-21)) {  # enough data for the sphere? # This is where it is failing
      meta_temp = as.matrix(meta_temp[-1,])
      # Select parts with no movement
      sdcriter = 0.013
      nomovement = which(meta_temp[,5] < sdcriter &
                           meta_temp[,6] < sdcriter &
                           meta_temp[,7] < sdcriter &
                           abs(as.numeric(meta_temp[,2])) < 2 &
                           abs(as.numeric(meta_temp[,3])) < 2 &
                           abs(as.numeric(meta_temp[,4])) < 2) # the latter three are to reduce chance of including clipping periods
      meta_temp = as.matrix(meta_temp[nomovement,])

      if (min(dim(meta_temp)) > 1) {
        meta_temp = as.matrix(meta_temp[(is.na(meta_temp[,4]) == F & is.na(meta_temp[,1]) == F),])
        npoints = nrow(meta_temp)
        cal.error.start = sqrt(as.numeric(meta_temp[,2])^2 + as.numeric(meta_temp[,3])^2 + as.numeric(meta_temp[,4])^2)
        cal.error.start = mean(abs(cal.error.start - 1))
        #check whether sphere is well populated
        tel = 0
        for (axis in 2:4) {
          if ( min(meta_temp[,axis]) < -spherecrit & max(meta_temp[,axis]) > spherecrit) {
            tel = tel + 1
          }
        }
        if (tel == 3) {
          spherepopulated = 1
        } else {
          spherepopulated = 0
          QC = "recalibration not done because not enough points on all sides of the sphere"
        }
      } else {
        cat("\nno non-movement found\n")
        QC = "recalibration not done because no non-movement data available"
        meta_temp = c()
      }
    } else {
      QC = "recalibration not done because not enough data in the file or because file is corrupt"
    }
    if (spherepopulated == 1) { #only try to improve calibration if there are enough datapoints around the sphere
      #---------------------------------------------------------------------------
      # START of Zhou Fang's code (slightly edited by vtv21 to use matrix meta_temp from above
      # instead the similar matrix generated by Zhou Fang's original code. This to allow for
      # more data to be used as meta_temp can now be based on 10 or more days of raw data
      input = meta_temp[,2:4] #as.matrix()
      if (use.temp == TRUE) { #at the moment i am always using temperature if mon == 2
        # mon == 2 & removed 19-11-2014 because it is redundant and would not allow for newer monitors to use it
        inputtemp = cbind(as.numeric(meta_temp[,8]),as.numeric(meta_temp[,8]),as.numeric(meta_temp[,8])) #temperature
      } else {
        inputtemp = matrix(0, nrow(input), ncol(input)) #temperature, here used as a dummy variable
      }
      meantemp = mean(as.numeric(inputtemp[,1]),na.rm=TRUE)
      inputtemp = inputtemp - meantemp
      offset = rep(0, ncol(input))
      scale = rep(1, ncol(input))
      tempoffset = rep(0, ncol(input))
      weights = rep(1, nrow(input))
      res = Inf
      maxiter = 1000
      tol = 1e-10
      for (iter in 1:maxiter) {
        curr = scale(input, center = -offset, scale = 1/scale) + scale(inputtemp, center = F, scale = 1/tempoffset)
        closestpoint = curr/ sqrt(rowSums(curr^2))
        k = 1
        offsetch = rep(0, ncol(input))
        scalech = rep(1,ncol(input))
        toffch = rep(0, ncol(inputtemp))
        for (k in 1:ncol(input)){
          fobj = lm.wfit(cbind(1, curr[,k],inputtemp[,k]) , closestpoint[,k, drop = F], w = weights)
          offsetch[k] = fobj$coef[1]
          scalech[k] = fobj$coef[2]
          if (use.temp == TRUE) {
            toffch[k] = fobj$coeff[3]
          }
          curr[,k] = fobj$fitted.values
        }
        offset = offset + offsetch / (scale  * scalech)
        if (use.temp == TRUE) {
          tempoffset = tempoffset * scalech + toffch
        }
        scale = scale * scalech
        res = c(res,  3 * mean(weights*(curr-closestpoint)^2/ sum(weights)))
        weights = pmin(1/ sqrt(rowSums((curr - closestpoint)^2)), 1/0.01)
        if (abs(res[iter+1] - res[iter]) < tol)  break
      }
      if (use.temp == FALSE) {
        meta_temp2 = scale(as.matrix(meta_temp[,2:4]),center = -offset, scale = 1/scale)
      } else {
        yy = as.matrix(cbind(as.numeric(meta_temp[,8]), as.numeric(meta_temp[,8]), as.numeric(meta_temp[,8])))
        meta_temp2 = scale(as.matrix(meta_temp[,2:4]),center = -offset, scale = 1/scale) +
          scale(yy, center = rep(meantemp,3), scale = 1/tempoffset)
      }     #equals: D2[,axis] = (D[,axis] + offset[axis]) / (1/scale[axis])
      # END of Zhou Fang's code
      #-------------------------------------------
      cal.error.end = sqrt(meta_temp2[,1]^2 + meta_temp2[,2]^2 + meta_temp2[,3]^2)
      cal.error.end = mean(abs(cal.error.end-1))
      # assess whether calibration error has sufficiently been improved
      if (cal.error.end < cal.error.start & cal.error.end < 0.01 & nhoursused > minloadcrit) { #do not change scaling if there is no evidence that calibration improves
        if (use.temp == TRUE) {
          QC = "recalibration done, no problems detected"
        } else if (use.temp == FALSE)  {
          QC = "recalibration done, but temperature values not used"
        }
        LD = 0 #stop loading
      } else { #continue loading data
        if (nhoursused > minloadcrit) {
          print(paste("new calibration error: ",cal.error.end," g",sep=""))
          print(paste("npoints around sphere: ", npoints,sep=""))
        }
        QC = "recalibration attempted with all available data, but possibly not good enough: Check calibration error variable to varify this"
      }
    }
    i = i + 1 #go to next block (12 hours-isch)
  } # End of While loop


  if (length(cal.error.end) > 0) {
    if (cal.error.end > cal.error.start) {
      QC = "recalibration not done because recalibration does not decrease error"
    }
  }
  if (length(ncol(meta_temp)) != 0) {
    spheredata = data.frame(A = meta_temp)
    if (use.temp == TRUE) {
      names(spheredata) = c("Euclidean Norm","meanx","meany","meanz","sdx","sdy","sdz","temperature")
    } else {
      names(spheredata) = c("Euclidean Norm","meanx","meany","meanz","sdx","sdy","sdz")
    }
  } else {
    spheredata = c()
  }
  QCmessage = QC
  if (printsummary == TRUE) {
    cat("\n----------------------------------------\n")
    cat("Summary of autocalibration procedure:\n")
    cat("\nStatus:\n")
    print(QCmessage)
    cat("\nCalibration error (g) before:\n")
    print(cal.error.start)
    cat("\nCallibration error (g) after:\n")
    print(cal.error.end)
    cat("\nOffset correction:\n")
    print(offset)
    cat("\nScale correction:\n")
    print(scale)
    cat("\nNumber of hours used:\n")
    print(nhoursused)
    cat("\nNumber of 10 second windows around the sphere:\n")
    print(npoints)
    cat("\nTemperature used (if available):\n")
    print(use.temp)
    cat("\nTemperature offset (if temperature is available):\n")
    print(tempoffset)
    cat("\n----------------------------------------\n")
  }
  if (use.temp==TRUE) {
    if (length(spheredata) > 0) {
      meantempcal = mean(spheredata[,8],na.rm=TRUE)
    } else {
      meantempcal = c()
    }
  } else {
    meantempcal = c()
  }
  invisible(list(scale=scale,offset=offset,tempoffset=tempoffset,
                 cal.error.start=cal.error.start,cal.error.end=cal.error.end,
                 spheredata=spheredata,npoints=npoints,nhoursused=nhoursused,
                 QCmessage = QCmessage,use.temp=use.temp,meantempcal=meantempcal,bsc_qc=bsc_qc))
}

#' @title g.downsample
#'
#' @name g.downsample
#'
#' @description Downsamples a vector of numeric values at three time resolutions:
#' 1 seconds, ws3 seconds, and ws2 second. Function is not intended for direct interaction
#' by package end user
#'
#' @param sig Vector of numeric values
#' @param fs Sample frequency
#' @param ws3 epoch size, e.g 5 seconds
#' @param ws2 epoch size, e.g 90 seconds
#'
#' @details List with three object: var1, var2, and var3 corresponding to downsample
#' time series at 1 seconds, ws2 seconds, and ws3 seconds resoluton, respectively
#'
#' @author Vincent T van Hees <vincentvanhees@gmail.com>
#'
#' @keywords internal

g.downsample = function(sig,fs,ws3,ws2) {
  #averaging per second => var1
  sig2 =cumsum(c(0,sig))
  select = seq(1,length(sig2),by=fs)
  var1 = diff(sig2[round(select)]) / abs(diff(round(select[1:(length(select))])))
  #averaging per ws3 => var2 (e.g. 5 seconds)
  select = seq(1,length(sig2),by=fs*ws3)
  var2 = diff(sig2[round(select)]) / abs(diff(round(select[1:(length(select))])))
  #averaging per ws2 => var3 (e.g. 15 minutes)
  select = seq(1,length(sig2),by=fs*ws2)
  var3 = diff(sig2[round(select)]) / abs(diff(round(select[1:(length(select))])))
  invisible(list(var1=var1,var2=var2,var3=var3))
}
