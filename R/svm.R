
#' @name svm
#'
#' @title sum vector magnitude
#'
#' @description  svm acts identically to 'mean.epoch', with the epoch set to the sampling period. In other words, it computes the
#' instantaneous sum of vector magnitudes of the acceleration at each record point. The function takes "AccData",
#' array and vector input. Note that if provided with an array with 4 or more columns, columns 2 to 4 are used -- the
#' first column is regard as a timestamp and hence ignored.
#'
#' @param obj AccData object
#' @param sqrt Function to use to calculate SVM
#'
#' @usage svm(obj, sqrt)
#'
#' @examples
#' dat <- read.bin(system.file("binfile/TESTfile.bin", package = "GENEAread")[1], calibrate = TRUE)
#' svm(dat)
#'
#' @export

svm <- function(obj,  sqrt){
  if (missing(sqrt)){
    sqrt = (length(dim(obj)) < 2)}
  if (length(dim(obj)) == 2) {
    obj = rowSums(obj[,-2:0 + min(ncol(obj), 4)]^2)
    if (sqrt) obj = sqrt(obj)
  } else {
    if (!sqrt) obj = obj^2
  }
  return(obj)
}
