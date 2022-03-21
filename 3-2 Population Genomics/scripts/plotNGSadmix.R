#!/bin/bash

## Plot NGSadmix results
## Adapted from https://github.com/ffertrindade/hz_toolkit for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-11-2022

# Reading arguments
args <- commandArgs(T)
if (length(args)!=3) {
  stop("Usage: Rscript ~/scripts/plotNGSadmix.R qopt_file k_num metadata", call.=FALSE)
}

qopt <- args[1]
k <- args[2]
metadata <- args[3]
rm(args)

# Get ID and pop info for each individual
name <- gsub(".qopt","",qopt)
popinfo <- read.csv(metadata, sep = ",", header = TRUE)
ind <- popinfo$run_accession
pop <- popinfo$Population
lat <- popinfo$latitude
lon <- popinfo$longitude

# Read inferred admixture proportions file
q <- read.table(qopt)

# Plot the results ordered by individual name
png(filename=paste("plots/",name, ".ordi.png", sep = ""), width=500, height=250)
ylab=paste("Admixture proportions (K=",k,")", sep = "")
main=paste(name,"(ordered by individuals name)", sep = " ")
par(mar=c(5,4,1,1))
barplot(t(q)[,order(ind)],
        col=c('darkred','blue4','darkgreen','orange','black','gray','pink'),
        space=0,
        border=NA,
        names=ind[order(ind)],
        las=2,
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
dev.off()

# Plot the results ordered by population name
png(filename=paste("plots/",name, ".ordp.png", sep = ""), width=500, height=500)
main=paste(name,"(ordered by populations name)", sep = " ")
par(mar=c(9,4,1,1))
barplot(t(q)[,order(pop)],
        col=c('darkred','blue4','darkgreen','orange','black','gray','pink'),
        space=0,
        border=NA,
        las=2,
        names=ind[order(pop)],
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
text(tapply(0.5:length(pop),pop[order(pop)],mean),-0.3,unique(pop[order(pop)]),xpd=T, srt=90, cex=0.75)
abline(v=cumsum(sapply(unique(pop[order(pop)]),function(x){sum(pop[order(pop)]==x)})),col=1,lwd=1)
dev.off()

# Plot the results ordered by latitude
png(filename=paste("plots/",name, ".lat.png", sep = ""), width=500, height=250)
main=paste(name,"(ordered by latitude)", sep = " ")
par(mar=c(5,4,1,1))
barplot(t(q)[,order(lat, decreasing = TRUE)],
        col=c('darkred','blue4','darkgreen','orange','black','gray','pink'),
        space=0,
        border=NA,
        las=2,
        names=ind[order(lat, decreasing = TRUE)],
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
dev.off()

# Plot the results ordered by longitude
png(filename=paste("plots/",name, ".lon.png", sep = ""), width=500, height=250)
main=paste(name,"(ordered by longitude)", sep = " ")
par(mar=c(5,4,1,1))
barplot(t(q)[,order(lon, decreasing = TRUE)],
        col=c('darkred','blue4','darkgreen','orange','black','gray','pink'),
        space=0,
        border=NA,
        las=2,
        names=ind[order(lon, decreasing = TRUE)],
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
dev.off()
