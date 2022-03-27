#!/bin/bash

## Construct the command line to run HaplotypeCaller for different samples, references and target regions
## Adapted from https://github.com/ffertrindade/proc_tools for the 3rd module Evolutionary Genomics of the curse Bioinformática y Genómica para la Biodiversidad
## Fernanda T. 03-07-2022

## keys and values
if [[ $# != 4 ]]; then
        echo "Usage: HaplotypeCaller.sh bamlist_file reference_file region dataset_name"
        exit 1
fi

BAMLIST=($(cat $1)) # List of bam files, ex. ~/leopard_data/leopard_17ind.bamlist
REFERENCE=$2 # Reference file, ex. ~/leopard_data/felcat9_chrE1.fasta
REGION=$3 # Region in the reference genome, ex. "NC_018736.3:1-10000000"
DATASET=$4 # Name of the dataset, ex. leopard_17ind

## creating command-line
REG=$(tr \: - <<<$REGION)
echo -ne "gatk3 \"-Xmx4g\" -T HaplotypeCaller -R $REFERENCE -L \"$REGION\" -ploidy 2 " > HaplotypeCaller.$REG.run.sh

for (( i=0; i<"${#BAMLIST[@]}"; i++ )); do
	echo -ne "-I ${BAMLIST[i]} " >> HaplotypeCaller.$REG.run.sh
done

echo -ne "> $DATASET.$REG.vcf" >> HaplotypeCaller.$REG.run.sh

bash HaplotypeCaller.$REG.run.sh 2> $DATASET.$REG.gatk.log

gzip $DATASET.$REG.vcf