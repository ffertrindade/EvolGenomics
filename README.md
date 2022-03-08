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
We'll use two different approaches to identify informative sites from the above mentioned dataset, which we'll call *leopard_17ind*: (i) genotype likelihoods and (ii) SNP calling. Note that depending on your dataset design (specially conserning depth) and the analyses and software you're planning to use, you'll choose one or another approach. In order to optimize run, you can check the [workflow](https://github.com/ffertrindade/EvolGenomics/main/day_8/workflows/calling), but here we're going to run step by step.

```
angsd -bam leopard_17ind.bamlist -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 24 -minQ 24 -SNP_pval 2e-6 -minMaf 0.05 -minInd 13 -out leopard_17ind -P 4
```
Below you can see how a beagle and maf file should look like:
![IMAGE HERE](http://url/to/img.png)

### Estimates of population structure and admixture

## Day 9: Adaptive evolution and demographic history
Today you are going to learn how to perform analyses of natural selection on complete genomes and how to infer historical demography. We're going to use as example the data from [Li et al., 2019](https://academic.oup.com/mbe/article/36/10/2111/5518928).

### Investigating adaptive evolution based on genome-wide datasets

