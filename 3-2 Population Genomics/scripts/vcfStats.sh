#!/bin/bash

## Run vcftools using a reduced VCF to get calling stats
## Adapted from https://speciationgenomics.github.io/filtering_vcfs/ for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-07-2022

## keys and values
if [[ $# != 2 ]]; then
	echo "Usage: vcfStats.sh input_vcf output_prefix"
        exit 1
fi

INPUT=$1 # VCF file, ex. leopard_17ind.filtered.vcf
OUT=$2 # Output prefix name, ex. leopard_17ind

### random sampling
bcftools view $INPUT | vcfrandomsample -r 0.012 > $OUT.subset.vcf
# compress vcf
bgzip $OUT.subset.vcf
# index vcf
bcftools index $OUT.subset.vcf.gz

### calculating stats

# mean depth per individual
vcftools --gzvcf $OUT.subset.vcf.gz --depth --out $OUT.subset
# proportion of missing data per individual
vcftools --gzvcf $OUT.subset.vcf.gz --missing-indv --out $OUT.subset
# homozygosity and the inbreeding coefficient F per individual
vcftools --gzvcf $OUT.subset.vcf.gz --het --out $OUT.subset

# mean depth per site
vcftools --gzvcf $OUT.subset.vcf.gz --site-mean-depth --out $OUT.subset
# proportion of missing data per site
vcftools --gzvcf $OUT.subset.vcf.gz --missing-site --out $OUT.subset
# quality score per site
vcftools --gzvcf $OUT.subset.vcf.gz --site-quality --out $OUT.subset
# allele frequency per variant
vcftools --gzvcf $OUT.subset.vcf.gz --freq2 --out $OUT.subset --max-alleles 2
