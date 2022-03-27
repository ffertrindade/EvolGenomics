#!/bin/bash

## Plot vcfStats results
## Adapted from https://speciationgenomics.github.io/filtering_vcfs/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-1-2022

## loading libraries and reading arguments
args <- commandArgs(T)
if (length(args)!=1) {
  stop("Usage: Rscript ~/scripts/vcfStats.R prefix_vcf_file", call.=FALSE)
}

input <- args[1]
rm(args)

library(tidyverse)

## Variant quality
filename <- paste(input, ".lqual", sep = "")
var_qual <- read_delim(filename, delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 500)
dev.off()

## Variant mean depth
filename <- paste(input, ".ldepth.mean", sep = "")
var_depth <- read_delim(filename, delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 10)
dev.off()

filename <- paste(input, ".ldepth.mean.sum.txt", sep = "")
sink(filename)
print(summary(var_depth$mean_depth))

filename <- paste(filename, "2.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 2)
dev.off()

filename <- paste(filename, "3.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_depth, aes(var_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 2)
dev.off()

## Variant missingness
filename <- paste(input, ".lmiss", sep = "")
var_miss <- read_delim(filename, delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

filename <- paste(input, ".lmiss.sum.txt", sep = "")
sink(filename)
print(summary(var_miss$fmiss))

filename <- paste(input, ".lmiss2.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

## Minor allele frequency
filename <- paste(input, ".frq", sep = "")
var_freq <- read_delim(filename, delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

filename <- paste(input, ".frq.sum.txt", sep = "")
sink(filename)
print(summary(var_freq$maf))

## Mean depth per individual
filename <- paste(input, ".idepth", sep = "")
ind_depth <- read_delim(filename, delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

## Proportion of missing data per individual
filename <- paste(input, ".imiss", sep = "")
ind_miss  <- read_delim(filename, delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

filename <- paste(filename, ".png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

## Heterozygosity and inbreeding coefficient per individual
filename <- paste(input, ".het", sep = "")
ind_het <- read_delim(filename, delim = "\t",
           col_names = c("ind","homObs", "homExp", "nsites", "f"), skip = 1)

filename <- paste(filename, ".f.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

filename <- paste(input, ".het.homObs.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(ind_het, aes(homObs/nsites)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + xlim(0, 1) +theme_light()
dev.off()

filename <- paste(input, ".het.homExp.png", sep = "")
png(filename = filename, width = 300, height = 300)
a <- ggplot(ind_het, aes(homExp/nsites)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + xlim(0, 1) + theme_light()
dev.off()
