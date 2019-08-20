
#' @name recalibrate
#'
#' @title recalibrate
#'
#' @description Taking a GENEActiv binfile and using the recalibration script to create a new calibrated binfile
#'
#' @param datadir The location of the directory/file containing GENEActiv binfile.
#' @param outputdir The location of the directory/file for the calibrated files to be saved.
#' @param use.temp Use temperature sensor data if available (Geneactive only)
#' @param spherecrit the minimum required acceleration value (in g) on both sides of 0 g for each
#' axis. Used to judge whether the sphere is sufficiently populated
#' @param minloadcrit the minimum number of hours the code needs to read for the autocalibration
#' procedure to be effective (only sensitive to multitudes of 12 hrs, other values will be ceiled).
#' After loading these hours only extra data is loaded if calibration
#' error has not been reduced to under 0.01 g.
#' @param printsummary if TRUE will print a summary when done chunksize number between 0.2 and 1
#' to specificy the size of chunks to be loaded as a fraction of a 12 hour period, e.g. 0.5 equals 6 hour chunks.
#'  The default is 1 (12 hrs). For machines with less than 4Gb of RAM memory a value below 1 is recommended.
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
#' @return Saves a calibrated binfile to an output folder
#'
#' @details Takes each binfile found in the data directory, calibrates according to the routine by Vincent T. van Hees
#' and saves the calibrated file to the specificied output directory
#'
#' @export
#'
#' @examples
#' \dontrun{
#' DataDirectory = "C:/Users/DataDirectory"
#' ReCalibrate(DataDirectory)
#' }

recalibrate = function(datadir,
                       outputdir,
                       use.temp = TRUE,
                       spherecrit = 0.3,
                       minloadcrit = 72,
                       printsummary = TRUE,
                       chunksize = c(0.5),
                       windowsizes = c(60,900,3600)){

    # Check that the datadir exists
    if (length(datadir) == 0) {
      if (length(datadir) == 0) {
        stop("Variable datadir is not defined")
      }
      if (missing(outputdir)) {
        print("Variable outputdir is not specified, creating an output directory.")
      }
    }

    filelist = FALSE

    # List all the files
    fnames=list.files(path = datadir)
    if (length(fnames) == 0){stop("There are no files in the data directory")}

    #### list of all bin files ####
    if (filelist == FALSE) {
      fnames = c(dir(datadir, recursive = TRUE, pattern="[.]bin"))
    }

    else {
      fnames = datadir
    }

    # If the outputdirectory is not specified create a defaulft
    if (missing(outputdir)) {
      outputfolder = paste(datadir,".Calibrated",sep="")
      dir.create(file.path(outputfolder))
    }

    for (i in 1:length(fnames)){

      Binfile = file.path(paste0(datadir, fnames[i]))

      # Find the calibration values of the data.
      C = GENEActiv.calibrate(Binfile,
                              use.temp = use.temp,
                              spherecrit = spherecrit,
                              minloadcrit = minloadcrit,
                              printsummary = printsummary,
                              chunksize = chunksize,
                              windowsizes = windowsizes)

      # Read in the current values from the file.
      Lines=readLines(Binfile,-1) # Reads the bin file to the point where the calibration data is.
      XOffset=Lines[49];  XGain=Lines[48]
      YOffset=Lines[51];  YGain=Lines[50]
      ZOffset=Lines[53];  ZGain=Lines[52]

      # Extracting just the numerical value. - Offset
      XOffsetN=as.numeric(unlist(strsplit(XOffset,"x offset:"))[2])
      YOffsetN=as.numeric(unlist(strsplit(YOffset,"y offset:"))[2])
      ZOffsetN=as.numeric(unlist(strsplit(ZOffset,"z offset:"))[2])

      # Extracting just the numerical value for the Gain
      XGainN=as.numeric(unlist(strsplit(XGain,"x gain:"))[2])
      YGainN=as.numeric(unlist(strsplit(YGain,"y gain:"))[2])
      ZGainN=as.numeric(unlist(strsplit(ZGain,"z gain:"))[2])

      # Calculating the New values using the calibration values- Rounding to whole numbers
      if (C$offset[1] != 0){
        XOffsetNew=round((XOffsetN * C$offset[1]) * 25.6)
        Lines[49]=paste("x offset:",sep="", XOffsetNew)}

      if (C$offset[2] != 0){
        YOffsetNew=round((YOffsetN * C$offset[2]) * 25.6)
        Lines[51]=paste("y offset:",sep="", YOffsetNew)
      }

      if (C$offset[2] != 0){
        ZOffsetNew=round((ZOffsetN * C$offset[3]) * 25.6)
        Lines[53]=paste("z offset:",sep="", ZOffsetNew)
      }

      if (C$scale[1] != 0){
        XGainNew=round(XGainN / C$scale[1])
        Lines[48]=paste("x gain:",sep="", XGainNew)
      }

      if (C$scale[2] != 0){
        YGainNew=round(YGainN / C$scale[2])
        Lines[50]=paste("y gain:",sep="", YGainNew)
        }
      if (C$scale[3] != 0){
        ZGainNew=round(ZGainN / C$scale[3])
        Lines[52]=paste("z gain:",sep="", ZGainNew)
      }

      # Creating the correct file name
      filename = fnames[i]
      names = strsplit(filename,".bin")
      writeLines(Lines, file.path(outputdir, paste0(names,"_Recalibrate.bin")))
    }
}

