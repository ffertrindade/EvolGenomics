#!/bin/bash

## Plot NGSadmix results
## Adapted from https://github.com/ffertrindade/hz_toolkit for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-11-2022

setwd("C:/Users/ffert/OneDrive - PUCRS - BR/Doutorado/analises/ngsadmix/")
setwd("G:/OneDrive - PUCRS - BR/Doutorado/analises/ngsadmix/")

# Get ID and pop info for each individual
name <- "batch001-005_198ind"
popinfo <- read.csv(paste(name, "_k2_sd1.popinfo", sep=""), sep = ";", header = FALSE)
ind <- popinfo$V1
pop <- popinfo$V2
reg <- popinfo$V3
lat <- popinfo$V4
lon <- popinfo$V5

# Read inferred admixture proportions file
q <- read.table(paste(name, "_k7_sd1.qopt", sep=""))

# Plot them (ordered by population or individual name)
ordq <- order(q[1])
ordr <- order(reg)

png(filename=paste("plots/",name, "_k7_sd1.ordi.png", sep = ""), width=2500, height=400)
par(mar=c(5,4,1,1))
barplot(t(q)[,order(ind)],
        col=c('darkred','blue4','darkgreen','black','gray','orange','pink'), #c(4,5,8,1,3)
        names=ind[order(ind)],
        las=2,
        ylab="Admixture proportions (K=7, seed=1)",
        main=name,
        cex.axis=1.5, 
        cex.names=0.7)
dev.off()

# order according to population only
png(filename=paste("plots/",name, "_k7_sd1.ordp.png", sep = ""), width=2500, height=400)
par(mar=c(8,4,1,1))
barplot(t(q)[,order(pop)],
        col=c('darkred','blue4','darkgreen','black','gray','orange','pink'), #c(4,5,8,1,3)
        space=0,
        border=NA,
        las=2,
        names=ind[order(pop)],
        ylab="Admixture proportions (K=7, seed=1)",
        main=name,
        cex.axis=1, 
        cex.names=0.7)
text(tapply(0.5:length(pop),pop[order(pop)],mean),-0.3,unique(pop[order(pop)]),xpd=T, srt=90, cex=0.75)
abline(v=cumsum(sapply(unique(pop[order(pop)]),function(x){sum(pop[order(pop)]==x)})),col=1,lwd=1)
dev.off()

# order according to lat only
png(filename=paste("plots/",name, "_k7_sd1.lat.png", sep = ""), width=2500, height=400)
par(mar=c(8,4,1,1))
barplot(t(q)[,order(lat, decreasing = TRUE)],
        col=c('darkred','blue4','darkgreen','black','gray','orange','pink'), #c(4,5,8,1,3)
        space=0,
        border=NA,
        las=2,
        names=ind[order(lat, decreasing = TRUE)],
        ylab="Admixture proportions (K=7, seed=1)",
        main=name,
        cex.axis=1, 
        cex.names=0.7)
dev.off()

# order according to lon only
png(filename=paste("plots/",name, "_k7_sd1.lon.png", sep = ""), width=2500, height=400)
par(mar=c(8,4,1,1))
barplot(t(q)[,order(lon, decreasing = TRUE)],
        col=c('darkred','blue4','darkgreen','black','gray','orange','pink'), #c(4,5,8,1,3)
        space=0,
        border=NA,
        las=2,
        names=ind[order(lon, decreasing = TRUE)],
        ylab="Admixture proportions (K=7, seed=1)",
        main=name,
        cex.axis=1, 
        cex.names=0.7)
dev.off()
