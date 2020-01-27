#!/bin/Rscript

stdError <- function(arg1,arg2) {
  return(sqrt((arg1-arg2)^2))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Input data file name not provided\n", call. = FALSE)
}

in_original <- args[1]
in_acr <- args[2]

heatTableOriginal <- read.table(in_original, col.names=c("xpos","ypos","temp"))
heatTableAcr <- read.table(in_acr, col.names=c("xpos","ypos","temp"))

if (!all.equal(heatTableOriginal[,c("xpos","ypos")],heatTableAcr[,c("xpos","ypos")])) {
  stop("The data grid of the two files does not match\n", call. = FALSE)
}

maxdens <- max(heatTableOriginal$temp)
mindens <- min(heatTableOriginal$temp)
range  <- stdError(maxdens,mindens)

summaryData <- data.frame(unclass(summary(stdError(heatTableAcr$temp, heatTableOriginal$temp) / (range))))
names(summaryData) <- ""
print(summaryData, digits = 3, zero.print = "0")
cat(sprintf("Range %.3f\n", range))
