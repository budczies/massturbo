IsMzWithinMargin <- function(value, mzValues, mzMaldiAdjustement = 1, mzMargin = 0.2) {
#  value <- myFilteredData$Mass
#  value <- tmpValuesThatShouldWork[1]
#  mzMargin = 0.2
#  mzMaldiAdjustement = 1

#  mzValuesAdjusted <- trunc(mzValues - mzMaldiAdjustement) # truncating is too strict..
  mzValuesAdjusted <- (mzValues - mzMaldiAdjustement)
  mzValuesAdjustedLowerEnd <- mzValuesAdjusted - mzMargin
  mzValuesAdjustedUpperEnd <- mzValuesAdjusted + mzMargin

#  isWithin <- FALSE
  isWithin <- "No"
  originalmzValue <- ""
  for (i in seq_along(mzValues)) {
    if ((value < mzValuesAdjustedUpperEnd[i]) & (value > mzValuesAdjustedLowerEnd[i])) {
#      isWithin <- TRUE
      isWithin <- "Yes"
      originalmzValue <- as.character(mzValues[i])
      break
    }
  }
#  return(isWithin)
  return(c(isWithin, originalmzValue))
}


#tmpValuesThatShouldWork <- myFilteredData[myFilteredData$Mass > 864 & myFilteredData$Mass < 865, -c(4:8) ]$Mass

#lapply(tmpValuesThatShouldWork, IsMzWithinMargin,  mzValues = mzValues, mzMaldiAdjustement = 1, mzMargin = 0.2)

