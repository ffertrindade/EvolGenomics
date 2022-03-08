#!/bin/bash

## Construct the command line to HaplotypeCaller for several samples and several target regions
## Adapted from https://github.com/ffertrindade/proc_tools for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-07-2022

## keys and values
if [[ $# != 4 ]]; then
        echo "Usage: HaplotypeCaller.sh bamlist_file reference_file region dataset_name"
        exit 1
fi

BAMLIST=($(cat $1)) # List of bam files, ex. leopard_17ind.bamlist
REFERENCE=$2 # Reference file, ex. felcat9_chrE1.fasta.gz
REGION=$3 # Region in the reference genome, ex. NC_018736.3:0-20000000
DATASET=$4 # Name of the dataset, ex. leopard_17ind

## creating command-line
REG=$(tr \: - <<<$REGION)
echo -e "gatk --java-options \"-Xmx4g\" HaplotypeCaller -R $REFERENCE -L \"$REGION\" -ploidy 2 -O $DATASET.$REGION.vcf.gz \ " > HaplotypeCaller.$REG.run.sh

for (( i=0; i<"${#BAMLIST[@]}"; i++ )); do
	echo -e "-I ${BAMLIST[i]} \ " >> HaplotypeCaller.$REG.run.sh
done

./HaplotypeCaller.$REG.run.sh