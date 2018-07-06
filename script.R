library(stringr)
library(openxlsx)

# which version
# which packages used
# as a command line
# which parameters accepts LCMS , MAldi, biomarker discovery, which need to comply with a format
# what output, table excel format with list of matches,

setwd("//Users//tanovsky//wip/Remy//HSIL")

# some parameter to be passed.. JUST DEMO VALUES
#mzValues <- c(1128.6754, 1365.7180,  915.4675, 1003.5491, 1032.5564, 1855.9615, 1875.9620)

# Remy check list
mzValues <- c(865.391, 901.556, 987.573, 1032.598, 1163.612, 1325.742, 1554.72, 1905.95)

# Just to test small subset..
#mzValues <- c(865.391, 2064.0052, 1809.8057, 1224.5485)

#mzValues <- c(865.391)

condition <- "HSIL"

#grepl(paste0("Intensity ", condition, "[0-9]+"), c("Intensity HSIL1", "sdf", "fdsfsd"))




#mzValues <- c(865.391, 865.391)

# Data file import
resourcesDir <- file.path(getwd(), "resources")
inputFilepath <- file.path(resourcesDir, "peptides.txt")
inputMembersHSIL_Filepath <- file.path(resourcesDir, "members HSIL.xlsx")

outputDir <- file.path(getwd(), "output")
outputFilepath <- file.path(outputDir, "resultTable.xlsx")
myData <- read.csv(file = inputFilepath, check.names = FALSE, sep = "\t")

membersHSILrawData <- read.xlsx(inputMembersHSIL_Filepath)
geneNameHSIL <- membersHSILrawData$Gene.name

columnNames <- colnames(myData)
pattern = "Intensity \\w"

conditionColumnAndCondition <- str_match(columnNames[grepl(pattern, columnNames)], "Intensity (?<condition>[a-z,A-Z]*)")
conditionColumn <- names(table(conditionColumnAndCondition[,1]))
numberOfSamplesPerCondition <- table(conditionColumnAndCondition[,2])
conditions <- names(numberOfSamplesPerCondition)


conditionColumnsHeaders <- columnNames[grepl(paste0("Intensity ", condition, "[0-9]+"), columnNames)]

selectedTableHeaders <- c("Mass", "Proteins", "Gene names", "Score", "Unique (Proteins)"
#                          , conditionColumnsHeaders
                          )

myFilteredData <- myData[, selectedTableHeaders]
#myFilteredData <- cbind(myFilteredData, paste0("#Positive", condition,"Replicates" = apply(myData[, conditionColumnsHeaders], 1, function(x) {sum(x > 0)})))
#paste0("#Positive", condition,"Replicates")







for (conditionColumnHeader in conditionColumn) {
  myFilteredData <- cbind(myFilteredData, apply(myData[, grepl(conditionColumnHeader, columnNames)], 1, mean))
}
colnames(myFilteredData) <- c(selectedTableHeaders, paste("Average", conditionColumn))


positiveReplicatesTable <- NULL
for (condition in conditions) {
  conditionColumnsHeaders <- columnNames[grepl(paste0("Intensity ", condition, "[0-9]+"), columnNames)]
#  tmpConditionColumnsHeaders <- columnNames[grepl(paste0("Intensity ", condition, "[0-9]+"), columnNames)]
  positiveReplicatesTable <- cbind(positiveReplicatesTable,  "#PositiveReplicates" = apply(myData[, conditionColumnsHeaders], 1, function(x) {sum(x > 0)}))
}
colnames(positiveReplicatesTable) <- paste0("#PositiveReplicates_", conditions)


myFilteredData <- cbind(myFilteredData, positiveReplicatesTable)

#isWithinRange <- sapply(myFilteredData$Mass, IsMzWithinMargin,  mzValues = mzValues, mzMaldiAdjustement = 1, mzMargin = 2)
#isWithinRange <- lapply(myFilteredData$Mass, IsMzWithinMargin,  mzValues = mzValues, mzMaldiAdjustement = 1, mzMargin = 0.4)
isWithinRange <- lapply(myFilteredData$Mass, IsMzWithinMargin,  mzValues = mzValues, mzMaldiAdjustement = 1, mzMargin = 0.2)

isWithinRange <- matrix(unlist(isWithinRange), ncol = 2, byrow = TRUE)
myFilteredData <- cbind(myFilteredData,
                        "mzWithinRange" = isWithinRange[, 1],
                        "MZ original" = isWithinRange[, 2],
                        "MaldiMass" = as.numeric(isWithinRange[, 2]) - 1,
                        "MaldiDelta" = (myFilteredData$Mass - (as.numeric(isWithinRange[, 2]) - 1)))



truncatedMasses <- trunc(myFilteredData$Mass)
uniqueTruncatedMasses <- sort(unique(truncatedMasses))

myFilteredData <- cbind(tmpTruncatedMass = truncatedMasses, myFilteredData)

filteredAndSortedTable <- NULL
for (i in uniqueTruncatedMasses) {
  tmpTable <- myFilteredData[myFilteredData$tmpTruncatedMass == i, -c(1)]
  filteredAndSortedTable <- rbind(filteredAndSortedTable, tmpTable[order(tmpTable$"Average Intensity HSIL"),])
}

filteredAndSortedTable <- cbind(filteredAndSortedTable,
                                genePresentAtHSIL = filteredAndSortedTable$"Gene names" %in% geneNameHSIL)


# To visualize only the interesting results..
filteredAndSortedTable[(filteredAndSortedTable$mzWithinRange != "No") &
                                                   filteredAndSortedTable$genePresentAtHSIL,]

dir.create(outputDir, showWarnings = FALSE)
write.xlsx(filteredAndSortedTable, file = outputFilepath)


# Microproteomic biomarker discovery-guided identification of MALDI imaging proteolytic peptides - a proof of concept illustrated with high grade dysplasia of the cervix.


