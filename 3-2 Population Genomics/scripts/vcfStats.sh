#!/bin/bash

## Run vcftools using a reduced VCF to get calling stats
## Adapted from https://github.com/ffertrindade/proc_tools for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-07-2022

## keys and values
if [[ $# != 2 ]]; then
	echo "Usage: vcfStats.sh input_vcf output_prefix"
        exit 1
fi

INPUT=$1 # VCF Input file, ex. ../allSamples_134ind_batch004.vcf.gz
OUT=$2 # Output prefix name, ex. allSamples_134ind_batch004

### random sampling
bcftools view $INPUT | vcfrandomsample -r 0.012 > $OUT.subset.vcf
# compress vcf
bgzip $OUT.subset.vcf
# index vcf
bcftools index $OUT.subset.vcf.gz

### calculating stats
vcftools --gzvcf $OUT.subset.vcf.gz --freq2 --out $OUT.subset --max-alleles 2
vcftools --gzvcf $OUT.subset.vcf.gz --depth --out $OUT.subset
vcftools --gzvcf $OUT.subset.vcf.gz --site-mean-depth --out $OUT.subset
vcftools --gzvcf $OUT.subset.vcf.gz --site-quality --out $OUT.subset
vcftools --gzvcf $OUT.subset.vcf.gz --missing-indv --out $OUT.subset
vcftools --gzvcf $OUT.subset.vcf.gz --missing-site --out $OUT.subset
vcftools --gzvcf $OUT.subset.vcf.gz --het --out $OUT.subset

