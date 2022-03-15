## 3-2 Population Genomics
Today you are going to learn how to perform basic populational analyses using complete genomes at low depth, including calling, diversity and homozygosity estimates, population structure and admixture. We're going to use as example the data from [Pecnerova et al., 2021](https://doi.org/10.1016/j.cub.2021.01.064).

The data we'll be using are **17 leopard individuals from 6 populations across Africa** (see [metadata](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind_metadata.csv)). These data had been previouslly subsampled and mapped to domestic cat [felCat9](https://www.ncbi.nlm.nih.gov/assembly/GCF_000181335.3) **chromosome E1** (63,494,689bp), and the resulted [bam](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind.bamlist) files are which we'll be using for this session. Bellow you can see basic stats of mapping quality and depth distribution for each individual along the chrE1. You can check the detailed workflow [here](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/workflows/mapping/snakefile).

### SNP calling vs. genotype likelihoods
We'll use two different approaches to identify informative sites from the above mentioned dataset, which we'll call *leopard_17ind*: (i) genotype likelihoods and (ii) SNP calling. Note that depending on your dataset design (specially conserning depth) and the analyses and software you're planning to use, you'll choose one or another approach. In order to optimize our run, you can check how we've created these files in [workflow](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/workflows/calling/snakefile), but here are the step by step in case you want to reproduce them (we'll give the files for the final exercises since the below steps can take too long).

The bam files are avaiable in a shared folder on the server. Please download the scripts from our repository:
```
mkdir scripts
cd scripts
curl -o HaplotypeCaller.sh https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/HaplotypeCaller.sh
curl -o vcfStats.sh https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/scripts/vcfStats.sh
curl -o plotNGSadmix.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotNGSadmix.R
cd ~
```

Estimate the genotype likelihoods using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD):
```
mkdir gl
cd gl
angsd -bam ~/leopard_data/leopard_17ind.bamlist -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 24 -minQ 24 -SNP_pval 2e-6 -minMaf 0.05 -minInd 17 -out leopard_17ind -P 4
cd ~
```
The above command resulted in 560,520 variable sites filtered according the parameters. Below you can see how our BEAGLE (file with genotype likelihoods for each individual) and MAF (file with allele frequencies for the entire population) files look like:
![BEAGLE](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_beagle.PNG)
![MAF](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_maf.PNG)

Perform SNP calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)
```
mkdir gatk
cd gatk
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:1-10000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:10000000-20000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:20000000-30000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:30000000-40000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:40000000-50000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:50000000-63494689" leopard_17ind
```
The above command resulted in a total of 1,448,909 identified sites. Bellow you can see how these VCF files looks like:
![GATK](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.NC_018736.3-20000000-30000000.vcf.gz.PNG)

However, the GATK command above results are too raw and need to be filtered. For that, we are going to use [vcftools](http://vcftools.sourceforge.net/man_latest.html), but there are a lot of options (even using GATK) to improve your dataset. Further, the bellow commands are very basic and general, to properlly filter your VCF files you can check them more deeply.
```
vcftools --gzvcf leopard_17ind.NC_018736.3-1-10000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-1-10000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --vcf leopard_17ind.NC_018736.3-10000000-20000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-10000000-20000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-20000000-30000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-20000000-30000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-30000000-40000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-30000000-40000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-40000000-50000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-40000000-50000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-50000000-63494689.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-50000000-63494689 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
bcftools concat -o leopard_17ind.filtered.vcf leopard_17ind.*.recode.vcf
```
The above command resulted in a total of 87,613 filtered sites. Bellow you can see how the final VCF file looks like:
![BCFTOOLS](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.filtered.vcf.PNG)

Futher, we can also use vcftools to retrieve some stats about our VCF file:
```
bash ~/scripts/vcfStats.sh leopard_17ind.filtered.vcf leopard_17ind
cd ~
```
Bellow are some of these stats:
![VCFSTATS]()

### Estimates of population structure and admixture
Now that we have SNPs and genotype likelihoods, the main exercise of today is explore genomic structure within our dataset. We're going to observe this by using [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix), which will give admixture proportions for each individual, and [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd), which will give (among other possible things) a matrix of clustering based in our dataset. Both of them use genotype likelihoods as input. 

Bellow is the command line to run NGSadmix for 2 ancestral populations. The analyses will create three files: (i) *leopard_17ind_k2.qopt*, which contains the ancestry proportions and should be used to plot the results; (ii) *leopard_17ind_k2.fopt.gz*, which contains calculated frequencies used to estimate the admixture proportions; (iii) *leopard_17ind_k2.log*, which contains usefull informations about the run (number of used sites, likelihood, etc). As an exercise, you should run it for K 2-6, compare the results and observe which K better describe our data based in the likelihood.

First, in case you haven't run the above commands to create beagle and VCF files, please download them from our repo:
```
mkdir gl
cd gl
curl -o leopard_17ind.beagle.gz https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/leopard_17ind.beagle.gz
cd ~
mkdir gatk
cd gatk
curl -o leopard_17ind.filtered.vcf.gz https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.filtered.vcf.gz
cd ~
```

And then, in order to run NGSadmix using the above beagle file, do:
```
mkdir ngsadmix
cd ngsadmix
NGSadmix -likes ../gl/leopard_17ind.beagle.gz -K 2 -P 4 -seed 1 -o leopard_17ind_k2
# Use the same command line above to run for K up to 6
```

Now we'll use a R script to plot our results. This script will plot the results ordering by individuals name, populations, latitude and longitude (see bellow). Observing these plots can better help to understand the population structure. Remember to plot the results for each K and, also as an excercise, compare the results and observe which K better describe our data based in the likelihood and biological interpretation. 
```
mkdir plots
Rscript ~/scripts/plotNGSadmix.R leopard_17ind_k2.qopt 2 ~/leopard_data/data/leopard_population_metadata.csv 
# Use the same command line above to run for K up to 6
```
Bellow are how the results for K2 ordered by population should look like:
![pop](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/ngsadmix/plots/leopard_17ind_k2.ordp.png)
