#!/bin/bash

## Plot Admixture results
## Adapted from https://speciationgenomics.github.io/ADMIXTURE/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-20-2022

## loading libraries and reading arguments
args <- commandArgs(T)
if (length(args)!=3) {
  stop("Usage: Rscript ~/scripts/plotAdmixture.R prefix_q_file k_num metadata", call.=FALSE)
}

input <- args[1] # prefix of Q file name, ex. leopard_17ind.filtered.pruned
k <- args[2] # Number of K from input Q file, ex. 2
metadata <- args[3] # Metadata ordered acording eigenvec file, ~/leopard_data/leopard_17ind_metadata.csv
rm(args)

library(ggplot2)

## Get ID and pop info for each individual
popinfo <- read.csv(metadata, sep = ",", header = TRUE)
ind <- popinfo$run_accession
pop <- popinfo$Population

## Read inferred admixture proportions file
q <- read.table(paste(input, ".", k, ".Q", sep = ""))

## Plot them (ordered by population or individual name)
ordp <- order(pop)

# plotting according individuals
png(filename=paste("plots/", input, ".", k, ".ordi.png", sep = ""), width=500, height=250)
ylab=paste("Admixture proportions (K=",k,")", sep = "")
main=paste(input,"(ordered by individuals name)", sep = " ")
par(mar=c(5,4,1,1))
barplot(t(q)[,order(ind)],
        col=c('darkred','blue4','orange','darkgreen','gray','black'),
        names=ind[order(ind)],
        las=2,
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
dev.off()

## plotting according population
png(filename=paste("plots/", input, ".", k, ".ordp.png", sep = ""), width=500, height=500)
ylab=paste("Admixture proportions (K=",k,")", sep = "")
main=paste(input,"(ordered by individuals name)", sep = " ")
par(mar=c(9,4,1,1))
barplot(t(q)[,order(pop)],
        col=c('darkred','blue4','orange','darkgreen','gray','black'), #c(4,5,8,1,3)
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
