# EvolGenomics
Tutorial of the practical activities for 3rd module *Evolutionary Genomics*  of the curse **Bioinformática y Genómica para la Biodiversidad**.

Module instructors:
- **Dr. Eduardo Eizirik**, Pontifícia Universidade Católica do Rio Grande do Sul (PUCRS), Brasil
- **Dr. Jorge Ramírez**, Universidad Nacional Mayor de San Marcos, Perú
- **Dr. Henrique Vieira Figueiró**, George Mason University, USA
- **Msc. Fernanda Trindade**, Pontifícia Universidade Católica do Rio Grande do Sul (PUCRS), Brasil

## Day 7: Phylogenomics
Today you are going to learn how to use complete genomes to reconstruct the groups phylogeny, to evaluate genealogical discordances and explore the role of hybridization in group's history. We're going to use as example the data from [Li et al., 2019](https://academic.oup.com/mbe/article/36/10/2111/5518928).

### Construct genomic windows
### Phylogenetic inference based on whole-genome sequences

## Day 8: Population Genomics
Today you are going to learn how to perform basic populational analyses using complete genomes at low depth, including calling, diversity and homozygosity estimates, population structure and admixture. We're going to use as example the data from [Pecnerova et al., 2021](https://doi.org/10.1016/j.cub.2021.01.064).

The data we'll be using are **17 leopard individuals from 6 populations across Africa** (see [metadata](https://github.com/ffertrindade/EvolGenomics/main/day_8/data/leopard_population_metadata.csv)). These data had been previouslly subsampled and mapped to domestic cat [felCat9](https://www.ncbi.nlm.nih.gov/assembly/GCF_000181335.3) **chromosome E1**, and the resulted [bam](https://github.com/ffertrindade/EvolGenomics/main/day_8/data/bam) files are which we'll be using for this session. Bellow you can see basic stats of mapping quality and depth distribution for each individual along the chrE1. You can check the detailed workflow [here](https://github.com/ffertrindade/EvolGenomics/main/day_8/workflows/mapping).
![IMAGE HERE](http://url/to/img.png)

### SNP calling vs. genotype likelihoods
We'll use two different approaches to identify informative sites from the above mentioned dataset, which we'll call *leopard_17ind*: (i) genotype likelihoods and (ii) SNP calling. Note that depending on your dataset design (specially conserning depth) and the analyses and software you're planning to use, you'll choose one or another approach. In order to optimize our run, you can check how we've created these files in [workflow](https://github.com/ffertrindade/EvolGenomics/main/day_8/workflows/calling), but here are the step by step in case you want to reproduce them (we'll give the files for the final exercises since the below steps take too long).

The bam files are avaiable in a shared folder on the server. Download the scripts from our repo:
```
mkdir scripts
wget https://github.com/ffertrindade/EvolGenomics/main/day_8/scripts
```
Estimate the genotype likelihoods using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD):
```
mkdir gl
cd gl
angsd -bam ../leopard_data/leopard_17ind.bamlist -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 24 -minQ 24 -SNP_pval 2e-6 -minMaf 0.05 -minInd 13 -setMaxDepth 85 -out leopard_17ind -P 4
cd ~
```
Below you can see how a beagle and maf file look like:
![IMAGE HERE](http://url/to/img.png)
Perform SNP calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)
```
mkdir gatk
cd gatk
../scripts/HaplotypeCaller.sh ../leopard_data/leopard_17ind.bamlist ../leopard_data/felcat9_chrE1.fasta.gz "NC_018736.3:1-20000000" leopard_17ind
../scripts/HaplotypeCaller.sh ../leopard_data/leopard_17ind.bamlist ../leopard_data/felcat9_chrE1.fasta.gz "NC_018736.3:20000000-40000000" leopard_17ind
../scripts/HaplotypeCaller.sh ../leopard_data/leopard_17ind.bamlist ../leopard_data/felcat9_chrE1.fasta.gz "NC_018736.3:40000000-63494689" leopard_17ind
vcftools --vcf leopard_17ind.NC_018736.3-1-20000000.vcf --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-1-20000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --vcf leopard_17ind.NC_018736.3-20000000-40000000.vcf --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-20000000-40000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --vcf leopard_17ind.NC_018736.3-40000000-63494689.vcf --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-40000000-63494689 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
bcftools concat -o leopard_17ind.filtered.vcf leopard_17ind.*.filtered.vcf
../scripts/vcfStats.sh leopard_17ind.filtered.vcf leopard_17ind
```

### Estimates of population structure and admixture

## Day 9: Adaptive evolution and demographic history
Today you are going to learn how to perform analyses of natural selection on complete genomes and how to infer historical demography. We're going to use as example the data from [Li et al., 2019](https://academic.oup.com/mbe/article/36/10/2111/5518928).

### Investigating adaptive evolution based on genome-wide datasets

