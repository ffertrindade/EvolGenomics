#!/bin/bash

## Plot PCA from plink results
## Adapted from https://speciationgenomics.github.io/pca/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-20-2022

## loading libraries and reading arguments
args <- commandArgs(T)
if (length(args)!=2) {
  stop("Usage: Rscript ~/scripts/plotPCA.R prefix_plink_file metadata", call.=FALSE)
}

input <- args[1]
metadata <- args[2] # ordered acording eigenvec file
rm(args)

library(ggplot2)

## read in data
pca <- read.csv(paste(input,".eigenvec",sep=""), sep = " ", header = FALSE)
eigenval <- scan((paste(input,".eigenval",sep="")))
popinfo <- read.csv(metadata, sep = ",", header = TRUE)

## sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# set populations
pops <- popinfo$Population
pca <- data.frame(pca, pops)

## convert to percentage variance explained
pve <- data.frame(PC = 1:17, pve = eigenval/sum(eigenval)*100)

## plot pca
png(filename=paste("plots/",input, ".PC1-2.png", sep = ""), width=500, height=250)
par(mar=c(5,4,1,1))
b <- ggplot(pca, aes(PC1, PC2, col = pops)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "orange", "gray", "black"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
dev.off()