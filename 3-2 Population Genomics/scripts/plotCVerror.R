#!/bin/bash

## Plot CV error of Admixture analyses
## Adapted from https://speciationgenomics.github.io/ADMIXTURE/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-20-2022

## loading libraries and reading arguments
args <- commandArgs(T)
if (length(args)!=1) {
  stop("Usage: Rscript ~/scripts/plotCVerror.R cverror_file", call.=FALSE)
}

input <- args[1]
rm(args)

library(ggplot2)

## Read CV error file file
cv_error <- read.csv(input, sep = " ", header = FALSE)

## plot cv error
png(filename=paste("plots/", input, ".png", sep = ""), width=500, height=500)
a <- ggplot(cv_error, aes(V1, V2)) + geom_bar(stat = "identity")
a + ylab("CV error") + xlab("K values") + theme_light()
dev.off()