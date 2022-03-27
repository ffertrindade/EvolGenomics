#!/bin/bash

## Plot PCAngsd results
## Adapted from http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/angsd-activity-pca-mds/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-20-2022

## loading libraries and reading arguments
args <- commandArgs(T)
if (length(args)!=3) {
  stop("Usage: Rscript ~/scripts/plotPCAngsd.R prefix_cov_file metadata", call.=FALSE)
}

input <- args[1]
metadata <- args[2] # ordered acording bamlist used to estimate the genotype likelihoods
rm(args)

library(ggplot2)

## read in the names of each individual and matrix file
popinfo <- read.csv(metadata, sep = ",", header = TRUE)
ind <- popinfo$run_accession
pop <- popinfo$Population
m <- as.matrix(read.table(paste(input, ".cov", sep="")))

## get the eigen vectors and edit data frame
e <- eigen(m)
row.names(e$vectors) <- ind
pca <- data.frame(e$vectors)
names(pca)[1:ncol(pca)] <- paste0("PC", 1:(ncol(pca)))
pca <- data.frame(pca, pop)

## get the eigen values
pve <- data.frame(PC = 1:17, pve = e[["values"]]/sum(e[["values"]])*100)

## plot the two first eigen vectors
png(filename=paste("plots/",input, ".PC1-2.png", sep = ""), width=500, height=250)
par(mar=c(5,4,1,1))
b <- ggplot(pca, aes(PC1, PC2, col = pop)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "orange", "gray", "black"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
dev.off()
