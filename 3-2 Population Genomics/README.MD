## 3-2 Population Genomics
Today you are going to learn how to perform basic populational analyses using complete genomes at low depth, including calling, basic statistics, and population structure using global ancestry proportions and principal component analyses. We're going to use as example the data from [Pecnerova et al., 2021](https://doi.org/10.1016/j.cub.2021.01.064).

![paper](https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/34083935-fab8-4250-9e3c-5f85f2bc009e/fx1.jpg)

The data we'll be using are **17 leopard individuals from 6 populations across Africa** (see [metadata](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind_metadata.csv)). These data had been previouslly subsampled and mapped to domestic cat [felCat9](https://www.ncbi.nlm.nih.gov/assembly/GCF_000181335.3) **chromosome E1** (63,494,689bp), and the resulted [bam](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/data/leopard_17ind.bamlist) files are which we'll be using for this session. You can check the detailed used workflow for mapping [here](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/workflows/mapping/snakefile).

### Practice on NLHPC cluster
Don't forget to create a script from **NLHPC cluster** with this [tool](https://wiki.nlhpc.cl/Generador_Scripts) to run your analyses, using [sbatch](https://wiki.nlhpc.cl/SISTEMA_GESTOR_DE_RECURSOS#Ejecutando_trabajos) to submit your jobs! But first, in case conda have not been previously configured for your user, do the following steps to do so: 
```
# 1- Login in your NLHPC user and enter your password
ssh instructor01@leftraru.nlhpc.cl
# 2- Load Anaconda3/2020.02 module
ml Anaconda3/2020.02
# 3- Activate conda for your account
conda init bash
# 4- Just logout and login again
# 5- Load Anaconda3/2020.02 module again
ml Anaconda3/2020.02
# 6- Activate conda bio-hpc environment (where all used software iss installed)
conda activate bio-hpc
```

Then, bellow is an example of how the script (to submit to the cluster with you jobs/analyses) should look like:
```
#!/bin/bash
#---------------SLURM Script - NLHPC ----------------
#SBATCH -J ngsadmix
#SBATCH -p general
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2400
#SBATCH --mail-user=fernanda.trindade@enorbita.org
#SBATCH --mail-type=ALL
#SBATCH -o ngsadmix_%j.out
#SBATCH -e ngsadmix_%j.err

#-----------------Toolchain---------------------------
# ----------------Modulos----------------------------
ml  Anaconda3/2020.02  
# ----------------Comando--------------------------
conda activate bio-hpc
NGSadmix -likes ../gl/leopard_17ind.beagle.gz -K 2 -P 4 -seed 1 -o leopard_17ind_k2
```
Remember that we also have **R scripts** to plot the results. You can download the scripts and the input files to your computer and then run using Rstudio (don't forget to install the required libraries), or run it directly in the cluster following the instructions.

## SNP calling vs. genotype likelihoods
From the above mentioned dataset, which we'll call *leopard_17ind*, we'll use two different approaches to identify informative sites: (i) **genotype likelihoods** (which is different from imputing or calling genotypes!) and (ii) **SNP calling**. Note that depending on your dataset design (specially concerning depth) and the analyses and software you're planning to use, you'll choose one or another approach. Further, you can use both approaches together to improve each other and, for example, do imputation. However, for today's class we are going to show how to estimate genotype likelihoods alone and how to call SNPs from a multi-sample perspective. Bellow here are the step by step in case you want to reproduce them (we'll give the files for the final exercises since the below steps can take too long).

The bam files are avaiable in a shared folder on the server. For the whole tutorial, we recomend you to download all the scripts and required files from our repository:
```
# Create and enter the directory
mkdir scripts
cd scripts
# Download all scripts using curl command and go to home directory
curl -o HaplotypeCaller.sh https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/HaplotypeCaller.sh
curl -o vcfStats.sh https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/vcfStats.sh
curl -o vcfStats.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/vcfStats.R
curl -o plotNGSadmix.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotNGSadmix.R
curl -o plotAdmixture.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotAdmixture.R
curl -o plotCVerror.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotCVerror.R
curl -o plotPCA.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotPCA.R
curl -o plotPCAngsd.R https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/scripts/plotPCAngsd.R
cd ~
# Create and enter the directory
mkdir leopard_data
cd leopard_data
# Download the required files for calling and metadata using curl command, unzip the reference with gunzip command, command and go to home directory
curl -o leopard_17ind_metadata.csv https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/leopard_17ind_metadata.csv
curl -o leopard_17ind.bamlist https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/leopard_17ind.bamlist # Remember to edit the correct path for the bam files in the server
curl -o felcat9_chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.fasta.gz
curl -o felcat9_chrE1.fasta.fai https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.fasta.fai
curl -o felcat9_chrE1.dict https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-2%20Population%20Genomics/data/felcat9_chrE1.dict
gunzip felcat9_chrE1.fasta.gz
```

### Genotype likelihoods using ANGSD
First, let's estimate the genotype likelihoods using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD):
```
# Create and enter the directory
mkdir gl
cd gl
# Run ANGSD (parameters discussed during class) and go to home directory
angsd -bam ~/leopard_data/leopard_17ind.bamlist -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 24 -minQ 24 -SNP_pval 2e-6 -minMaf 0.05 -minInd 17 -out leopard_17ind -P 4
cd ~
```
**Parameters used:**
```
-GL [int]: Estimate genotype likelhoods from the mapped reads. Four different methods are available.
        1: SAMtools model
        2: GATK model
        3: SOAPsnp model
        4: SYK model
-doMajorMinor [int]: Infer Major and MInor alleles.
	    1: Infer major and minor from GL
	    2: Infer major and minor from allele counts
	    3: use major and minor from a file (requires -sites file.txt)
	    4: Use reference allele as major (requires -ref)
	    5: Use ancestral allele as major (requires -anc)
-doMaf [int]: Calculate persite frequencies '.mafs.gz'.
        1: Frequency (fixed major and minor). Here both the major and minor allele is assumed to be known (inferred or given by user). The allele frequency is the obtained using based on the genotype likelihoods. 
        2: Frequency (fixed major unknown minor). Here the major allele is assumed to be known (inferred or given by user) however the minor allele is not determined. Instead we sum over the 3 possible minor alleles weighted by their probabilities.
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
-doGlf [int]: Output the log genotype likelihoods to a file.
        1: binary all 10 log genotype likelihood
        2: beagle genotype likelihood format (use directly for imputation)
        3: beagle binary
        4: textoutput of all 10 log genotype likelihoods.
-minMapQ [int]: Minimum mapQ quality.
-minQ [int]: Minimum base quality score.
-SNP_pval [float]: Only work with sites with a p-value less than [float]. The p-value used for calling snps. SNPs are called based on their allele frequencies. If a site has a minor allele frequency significantly different from 0 a site is called as polymorphic.
-minMaf [float]: Only work with sites with a maf above [float].
-minInd [int]: Only keep sites with at least 1x from at least [int] individuals.
```
The above command resulted in **560,520 variable sites filtered according the parameters**. Below you can see how our BEAGLE (file with genotype likelihoods for each individual) and MAF (file with allele frequencies for the entire population) files look like:
![BEAGLE](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_beagle.PNG)
![MAF](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/gl/file_maf.PNG)

### SNP calling using GATK and VCFTOOLS
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
**Parameters used:**
```
-T HaplotypeCaller : method used.
-R : reference fasta file (must have .dict and.fai).
-L : search only in the region of interest.
-ploidy 2 : diploid.
-I : input bam file (several in the same run).
# Note that can be more useful to use GVCF approach https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
```
The above command resulted in a total of **1,448,909 identified sites**. 

However, the GATK command above results are too raw and need to be filtered. For that, we are going to use [vcftools](http://vcftools.sourceforge.net/man_latest.html), but there are a lot of options (even using GATK) to improve your dataset. Further, the bellow commands are very basic and general, to properlly filter your VCF files you can check them more deeply.
```
# Run vcftools (parameters discussed during class) to filter the raw VCFs
vcftools --gzvcf leopard_17ind.NC_018736.3-1-10000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-1-10000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-10000000-20000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-10000000-20000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-20000000-30000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-20000000-30000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-30000000-40000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-30000000-40000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-40000000-50000000.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-40000000-50000000 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
vcftools --gzvcf leopard_17ind.NC_018736.3-50000000-63494689.vcf.gz --recode --recode-INFO-all --out leopard_17ind.NC_018736.3-50000000-63494689 --minDP 5 --maxDP 15 --maf .05 --max-maf .95 --min-alleles 2 --max-alleles 2 --remove-indels
# Concatenate the filtered VCFs created with vcftools using bcftools
bcftools concat -o leopard_17ind.filtered.vcf leopard_17ind.NC_018736.3-1-10000000.recode.vcf leopard_17ind.NC_018736.3-10000000-20000000.recode.vcf leopard_17ind.NC_018736.3-20000000-30000000.recode.vcf leopard_17ind.NC_018736.3-30000000-40000000.recode.vcf leopard_17ind.NC_018736.3-40000000-50000000.recode.vcf leopard_17ind.NC_018736.3-50000000-63494689.recode.vcf
```
**Parameters used:**
```
--gzvcf : input gziped VCF file.
--out : output prefix name.
--recode and --recode-INFO-all : recode header information.
--minDP <float> and --maxDP <float> : Includes only genotypes greater than or equal to the "--minDP" value and less than or equal to the "--maxDP" value. This option requires that the "DP" FORMAT tag is specified for all sites.
--maf <float> and --max-maf <float> : Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value and less than or equal to the "--max-maf" value. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site.
--min-alleles <integer> and --max-alleles <integer> : Include only sites with a number of alleles greater than or equal to the "--min-alleles" value and less than or equal to the "--max-alleles" value. 
--remove-indels : Exclude sites that contain an indel. For these options "indel" means any variant that alters the length of the REF allele.
```
The above command resulted in a total of **105,321 filtered sites**. Bellow you can see how the final VCF file looks like:
![BCFTOOLS](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/leopard_17ind.filtered.vcf.png)

Finally, to have a stats summary of our cleanned VCF (which is very useful to understand our data and exclude outliers) we can also use vcftools:
```
# Run vcftools using the downloaded script (shown during class) and go to home directory
bash ~/scripts/vcfStats.sh leopard_17ind.filtered.vcf leopard_17ind
Rscript ~/scripts/vcfStats.R leopard_17ind.subset
cd ~
```

Bellow are some of these stats (see [here](https://github.com/ffertrindade/EvolGenomics/tree/main/3-2%20Population%20Genomics/results/vcf/vcfstats) for more):

![IDEPTH](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/vcfstats/leopard_17ind.subset.idepth.png)
![LQUAL](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/vcf/vcfstats/leopard_17ind.subset.lqual.png)

## Estimates of population structure and admixture
Now that we have SNPs and genotype likelihoods, the main exercise of today is explore genomic structure within our dataset. We're going to observe this by using [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) and [Admixture](https://dalexander.github.io/admixture/), which will give global admixture proportions for each individual, and [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) and [PLINK](https://zzz.bwh.harvard.edu/plink/), which will give (and can give other several things) a matrix of clustering based in our dataset. Those four software will be using the **beagle or VCF** files. 

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

### Global ancestry proportions using NGSadmix
Bellow is the command line to run **NGSadmix** for 2 ancestral populations based on beagle file. The analyses will create three files: (i) *leopard_17ind_k2.qopt*, which contains the ancestry proportions and should be used to plot the results; (ii) *leopard_17ind_k2.fopt.gz*, which contains calculated frequencies used to estimate the admixture proportions; (iii) *leopard_17ind_k2.log*, which contains usefull informations about the run (number of used sites, likelihood, etc). As an **exercise of the course, you should run it for K 2-6**, compare the results and observe which K better describe our data based in the likelihood.
```
# Create and enter the directory
mkdir ngsadmix
cd ngsadmix
# Run NGSadmix (parameters discussed during class) in order to get individuals ancestry proportions
NGSadmix -likes ../gl/leopard_17ind.beagle.gz -K 2 -P 4 -seed 1 -o leopard_17ind_k2
## Use the same command line above to run for K up to 6 ##
```
**Parameters used:**
```
-likes : input beagle GL file.
-K : number of ancestral populations to be modeled.
-P : number of threads to be used.
-seed : number seed to reproduce the run.
-o : output prefix name.
```

Now we'll use a R script to plot our results. This script will plot the results ordering by individuals name, populations, latitude and longitude (accordint to the metadata file). Observing these plots can better help to understand the population structure. **Remember to plot the results for each K** and, also as an excercise, compare the results and observe which K better describe our data based in the likelihood and biological interpretation. 
```
# Create the directory where the plots will be saved
mkdir plots
# Plot NGSadmix results using the downloaded script (shown during class) and go to home directory
Rscript ~/scripts/plotNGSadmix.R leopard_17ind_k2.qopt 2  ~/leopard_data/leopard_17ind_metadata.csv 
## Use the same command line above to plot for K up to 6 ##
cd ~
```

### PCA using PLINK
Now we're going to use the VCF file to observe population structure and clustering. First we need to filter linked sites using **PLINK**. The analyses will create four files: (i) *.prune.in*, list of unlinked sites; (ii) *.prune.out*, list of linked sites to exclude; (vi) *.log*, log file of PLINK run; (vii) *.nosex*, individuals with unknown sex. As an **exercise of the course, you should run it and do the PCA plot**, and compare the results with the other analyses.

Therefore, let's filtering our VCF to contain only unlinked sites:
```
# Enter the directory
cd gatk
# Run PLINK (parameters discussed during class) to get a list of unlinked sites and go to home directory
plink --vcf leopard_17ind.filtered.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 10 10 0.9 --out leopard_17ind.filtered
cd ~
```
**Parameters used:**
```
--vcf : input VCF file.
--double-id : when you dont have pedigree/family data, then it'll just double the IDs' names.
--allow-extra-chr : allow different chromosomes than expected in humans.
--set-missing-var-ids @:# : since our SNPs are not annotated it'll stands "chromosome:position" as SNP ID
--indep-pairwise 10 10 0.9 : windows of 10 SNPs, step size of 10 SNPs, prune sites with r2 greater than 0.9
--out : output prefix name.
```
The above command shows that we have **39,666 unlinked sites** in this dataset.

Bellow is the command line to run **PLINK** to perform PCA based on VCF file. The analyses will create seven files: (i) *.bed*, binary file with genotypes, it'll be used for next Admixture analysis; (ii) *.bim*, file with a sort of information about the sites; (iii) *.eigenval*, eigenvalues for the PCA analysis; (iv) *.eigenvec*, eigenvectors for the PCA analysis; (v) *.fam*, file with information about the individuals; (vi) *.log*, log file of PLINK run; (vii) *.nosex*, individuals with unknown sex. As an **exercise of the course, you should run it**, and compare the results with the other analyses:
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
**Parameters used:**
```
--extract : use only this list of sites.
--make-bed : create the bed, bim and fam files as output.
--pca : calculate eigenvalues and eigenvectors.
```

### Global ancestry proportions using Admixture
The bellow steps are not required to perform for the course exercises, but you can do if you'd like.

We can now use some of the other files created by PLINK to run **Admixture** and get the individual's global ancestry proportions. This analysis is compared to that from NGSadmix, but using different dataset and premisses. Do the following:
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
The final plots should look like [these](https://github.com/ffertrindade/EvolGenomics/tree/main/3-2%20Population%20Genomics/results/admixture/plots).

![ADMIXTURE](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/admixture/plots/leopard_17ind.filtered.pruned.2.ordp.png)

### PCA using PCAngsd
Finally, you can also run **PCAngsd** using your beagle file. Since it's a script, you should be able to download it from github:
```
# Get the repository from github (instalation according http://www.popgen.dk/software/index.php/PCAngsd) inside your scripts directory
cd scripts
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
python setup.py build_ext --inplace
pip install --user -r requirements.txt
cd ~
```

And then you can perform the analyses as follow:
```
# Create and enter the directory
mkdir pcangsd
cd pcangsd
# Run PCAngsd using you beagle file
python ~/scripts/pcangsd/pcangsd.py -beagle ../gl/leopard_17ind.beagle.gz -out leopard_17ind -threads 4
# Plot PCAngsd results using the previously downloaded script
Rscript ~/scripts/plotPCAngsd.R leopard_17ind ~/leopard_data/leopard_17ind_metadata.csv
```
Here is the expected result:
![PCANGSD](https://github.com/ffertrindade/EvolGenomics/blob/main/3-2%20Population%20Genomics/results/pcangsd/plots/leopard_17ind.PC1-2.png)