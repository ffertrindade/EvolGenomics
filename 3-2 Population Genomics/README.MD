## 3-2 Population Genomics
Today you are going to learn how to perform basic populational analyses using complete genomes at low depth, including calling, diversity and homozygosity estimates, population structure and admixture. We're going to use as example the data from [Pecnerova et al., 2021](https://doi.org/10.1016/j.cub.2021.01.064).

The data we'll be using are **17 leopard individuals from 6 populations across Africa** (see [metadata](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind_metadata.csv)). These data had been previouslly subsampled and mapped to domestic cat [felCat9](https://www.ncbi.nlm.nih.gov/assembly/GCF_000181335.3) **chromosome E1** (63,494,689bp), and the resulted [bam](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind.bamlist) files are which we'll be using for this session. Bellow you can see basic stats of mapping quality and depth distribution for each individual along the chrE1. You can check the detailed workflow [here](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/workflows/mapping/snakefile).

### SNP calling vs. genotype likelihoods
We'll use two different approaches to identify informative sites from the above mentioned dataset, which we'll call *leopard_17ind*: (i) genotype likelihoods and (ii) SNP calling. Note that depending on your dataset design (specially conserning depth) and the analyses and software you're planning to use, you'll choose one or another approach. In order to optimize our run, you can check how we've created these files in [workflow](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/workflows/calling/snakefile), but here are the step by step in case you want to reproduce them (we'll give the files for the final exercises since the below steps can take too long).

The bam files are avaiable in a shared folder on the server. Please download the scripts from our repository:
```
# Create and enter the directory
mkdir scripts
cd scripts
# Download all scripts using curl command and go to home directory
curl -o HaplotypeCaller.sh https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/HaplotypeCaller.sh
curl -o vcfStats.sh https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/scripts/vcfStats.sh
curl -o plotNGSadmix.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotNGSadmix.R
curl -o vcfStats.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/vcfStats.R
curl -o plotAdmixture.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotAdmixture.R
curl -o plotCVerror.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotCVerror.R
curl -o plotPCA.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotPCA.R
curl -o plotPCAngsd.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotPCAngsd.R
cd ~
# Create and enter the directory
mkdir leopard_data
cd leopard_data
# Download all scripts using curl, unzip the reference with gunzip command, command and go to home directory
curl -o leopard_17ind_metadata.csv https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/leopard_17ind_metadata.csv
curl -o leopard_17ind.bamlist https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/leopard_17ind.bamlist
curl -o felcat9_chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.fasta.gz
curl -o felcat9_chrE1.fasta.fai https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.fasta.fai
curl -o felcat9_chrE1.dict https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.dict
gunzip felcat9_chrE1.fasta.gz
```

Estimate the genotype likelihoods using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD):
```
# Create and enter the directory
mkdir gl
cd gl
# Run ANGSD (parameters discussed during class) and go to home directory
angsd -bam ~/leopard_data/leopard_17ind.bamlist -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 24 -minQ 24 -SNP_pval 2e-6 -minMaf 0.05 -minInd 17 -out leopard_17ind -P 4
cd ~
```
The above command resulted in **560,520 variable sites filtered according the parameters**. Below you can see how our BEAGLE (file with genotype likelihoods for each individual) and MAF (file with allele frequencies for the entire population) files look like:
![BEAGLE](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_beagle.PNG)
![MAF](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_maf.PNG)

Perform SNP calling (creating a VCF file) using [GATK](https://gatk.broadinstitute.org/hc/en-us):
```
# Create and enter the directory
mkdir gatk
cd gatk
# Run GATK using the downloaded script (shown during class) for each region of the E1 chromosome
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:1-10000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:10000000-20000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:20000000-30000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:30000000-40000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:40000000-50000000" leopard_17ind
bash ~/scripts/HaplotypeCaller.sh ~/leopard_data/leopard_17ind.bamlist ~/leopard_data/felcat9_chrE1.fasta "NC_018736.3:50000000-63494689" leopard_17ind
```
The above command resulted in a total of **1,448,909 identified sites**. Bellow you can see how these VCF files looks like:
![GATK](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.NC_018736.3-20000000-30000000.vcf.gz.PNG)

However, the GATK command above results are too raw and need to be filtered. For that, we are going to use [vcftools](http://vcftools.sourceforge.net/man_latest.html), but there are a lot of options (even using GATK) to improve your dataset. Further, the bellow commands are very basic and general, to properlly filter your VCF files you can check them more deeply.
```
# Run vcftools (parameters discussed during class) to filter the raw VCFs
vcftools --gzvcf leopard_17ind.NC_018736.3-1-10000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-1-10000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gvcf leopard_17ind.NC_018736.3-10000000-20000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-10000000-20000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-20000000-30000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-20000000-30000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-30000000-40000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-30000000-40000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-40000000-50000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-40000000-50000000 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-50000000-63494689.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-50000000-63494689 --minDP 5 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
# Concatenate the filtered VCFs created with vcftools using bcftools
bcftools concat -o leopard_17ind.filtered.vcf leopard_17ind.*.recode.vcf
```
The above command resulted in a total of **87,613 filtered sites**. Bellow you can see how the final VCF file looks like:
![BCFTOOLS](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.filtered.vcf.PNG)

Futher, to have a stats summary of our cleanned VCF (which is very useful to understand our data and exclude outliers) we can also use vcftools:
```
# Run vcftools using the downloaded script (shown during class) and go to home directory
bash ~/scripts/vcfStats.sh leopard_17ind.filtered.vcf leopard_17ind
Rscript ~/scripts/vcfStats.R leopard_17ind.subset
cd ~
```

Bellow are some of these stats (see [here](https://github.com/ffertrindade/EvolGenomics/tree/main/3-2%20Population%20Genomics/results/vcf/vcfstats) for more):

![IDEPTH](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/vcfstats/leopard_17ind.subset.idepth.png)
![LQUAL](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/vcfstats/leopard_17ind.subset.lqual.png)

### Estimates of population structure and admixture
Now that we have SNPs and genotype likelihoods, the main exercise of today is explore genomic structure within our dataset. We're going to observe this by using [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) and [Admixture](https://dalexander.github.io/admixture/), which will give admixture proportions for each individual, and [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) and [PLINK](https://zzz.bwh.harvard.edu/plink/), which will give (and can give other several things) a matrix of clustering based in our dataset. Those four software will be using the **beagle or VCF** files. 

First of all, **in case you haven't run the above commands to create beagle and VCF files, please download them** from our repo:
```
# Create and enter the directory
mkdir gl
cd gl
# Download the beagle file using curl command and go to home directory
curl -o leopard_17ind.beagle.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/results/gl/leopard_17ind.beagle.gz
cd ~
# Create and enter the directory
mkdir gatk
cd gatk
# Download the filtered VCF file using curl command, unzip it using gunzip command, and go to home directory
curl -o leopard_17ind.filtered.vcf.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.filtered.vcf.gz
gunzip leopard_17ind.filtered.vcf.gz
cd ~
```

Bellow is the command line to run **NGSadmix** for 2 ancestral populations based on beagle file. The analyses will create three files: (i) *leopard_17ind_k2.qopt*, which contains the ancestry proportions and should be used to plot the results; (ii) *leopard_17ind_k2.fopt.gz*, which contains calculated frequencies used to estimate the admixture proportions; (iii) *leopard_17ind_k2.log*, which contains usefull informations about the run (number of used sites, likelihood, etc). As an **exercise of the course, you should run it for K 2-6**, compare the results and observe which K better describe our data based in the likelihood.
```
# Create and enter the directory
mkdir ngsadmix
cd ngsadmix
# Run NGSadmix (parameters discussed during class) in order to get individuals ancestry proportions
NGSadmix -likes ../gl/leopard_17ind.beagle.gz -K 2 -P 4 -seed 1 -o leopard_17ind_k2
## Use the same command line above to run for K up to 6 ##
```

Now we'll use a R script to plot our results. This script will plot the results ordering by individuals name, populations, latitude and longitude (see bellow). Observing these plots can better help to understand the population structure. Remember to plot the results for each K and, also as an excercise, compare the results and observe which K better describe our data based in the likelihood and biological interpretation. 
```
# Create the directory where the plots will be saved
mkdir plots
# Plot NGSadmix results using the downloaded script (shown during class) and go to home directory
Rscript ~/scripts/plotNGSadmix.R leopard_17ind_k2.qopt 2  ~/leopard_data/leopard_17ind_metadata.csv 
## Use the same command line above to plot for K up to 6 ##
cd ~
```
Bellow are how the results for K2 ordered by population should look like:
![pop](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/ngsadmix/plots/leopard_17ind_k2.ordp.png)

Now we're going to use the VCF file to observe population structure and clustering. First we need to filter unlinked sites using **PLINK**. The analyses will create four files: (i) *.in*, list of unliked sites; (ii) *.out*, list of liked sites to exclude; (vi) *.log*, log file of PLINK run; (vii) *.nosex*, individuals with unknown sex. As an **exercise of the course, you should run it up to PCA plots**, and compare the results with the other analyses.

Therefore, lets filtering our VCF to contain only unlinked sites:
```
# Enter the directory
cd gatk
# Run PLINK (parameters discussed during class) to get a list of unlinked sites and go to home directory
plink --vcf leopard_17ind.filtered.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 10 10 0.9 --out leopard_17ind.filtered
cd ~
```
The above command shows that we have **34,497 unlinked sites** in this dataset.

Bellow is the command line to run **PLINK** to perform PCA based on VCF file. The analyses will create seven files: (i) *.bed*, binary file with genotypes, it'll be used for next Admixture analysis; (ii) *.bim*, file with information about the sites; (iii) *.eigenval*, eigenvalues for the PCA analysis; (iv) *.eigenvec*, eigenvectors for the PCA analysis; (v) *.fam*, file with information about the individuals; (vi) *.log*, log file of PLINK run; (vii) *.nosex*, individuals with unknown sex. As an **exercise of the course, you should run it**, and compare the results with the other analyses:
```
# Create and enter the directory
mkdir pca
cd pca
# Run PLINK (parameters discussed during class) to perform PCA analysis
plink --vcf ../gatk/leopard_17ind.filtered.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ../gatk/leopard_17ind.filtered.prune.in --make-bed --pca --out leopard_17ind.filtered.pruned
# Create the directory where the plots will be saved
mkdir plots
# Plot PCA results using the downloaded script (shown during class) and go to home directory
Rscript ~/scripts/plotPCA.R leopard_17ind.filtered.pruned ~/leopard_data/leopard_17ind_metadata.csv
cd ~
```

We can now use some of the files created by PLINK to run **Admixture** and get the individual's ancestry proportions. This analysis is capared to that from NGSadmix, but using different dataset and premisses. Do the following:
```
# Create and enter the directory
mkdir admixture
cd admixture
# Edit bim file to be properlly used by Admixture and copy the necessary files to the directory where you will run it
awk '{$1="0";print $0}' ../pca/leopard_17ind.filtered.pruned.bim > leopard_17ind.filtered.pruned.bim
cp ../pca/leopard_17ind.filtered.pruned.bed ../pca/leopard_17ind.filtered.pruned.fam .
# Run Admixture (parameters discussed during class) for several K values
admixture --cv leopard_17ind.filtered.pruned.bed 2 > leopard_17ind.filtered.pruned.k2.log
admixture --cv leopard_17ind.filtered.pruned.bed 3 > leopard_17ind.filtered.pruned.k3.log
admixture --cv leopard_17ind.filtered.pruned.bed 4 > leopard_17ind.filtered.pruned.k4.log
admixture --cv leopard_17ind.filtered.pruned.bed 5 > leopard_17ind.filtered.pruned.k5.log
admixture --cv leopard_17ind.filtered.pruned.bed 6 > leopard_17ind.filtered.pruned.k6.log
# retrieve CV error for all runs
grep "CV error" leopard_17ind.filtered.pruned.k*.log | awk '{print $3,$4}' | sed 's/[(|)|=|:]//g' > leopard_17ind.filtered.pruned.cverror.txt
# Create the directory where the plots will be saved
mkdir plots
# Plot the CV error distribution using the downloaded script
Rscript ~/scripts/plotCVerror.R leopard_17ind.filtered.pruned.cverror.txt
# Plot the Admixture results for all K using the downloaded script
Rscript ~/scripts/plotAdmixture.R leopard_17ind.filtered.pruned 2 ~/leopard_data/leopard_17ind_metadata.csv
Rscript ~/scripts/plotAdmixture.R leopard_17ind.filtered.pruned 3 ~/leopard_data/leopard_17ind_metadata.csv
Rscript ~/scripts/plotAdmixture.R leopard_17ind.filtered.pruned 4 ~/leopard_data/leopard_17ind_metadata.csv
Rscript ~/scripts/plotAdmixture.R leopard_17ind.filtered.pruned 5 ~/leopard_data/leopard_17ind_metadata.csv
Rscript ~/scripts/plotAdmixture.R leopard_17ind.filtered.pruned 6 ~/leopard_data/leopard_17ind_metadata.csv
cd ~
```