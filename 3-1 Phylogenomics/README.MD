## 3-1 Phylogenomics
Our dataset was obtained from our 2019 paper [Li et al., 2019](https://academic.oup.com/mbe/article/36/10/2111/5518928). For this tutorial, we have selected chromosome E1 (~61Mb) from six cat species: Lynx canadensis, Lynx rufus, Lynx lynx, Panthera pardus,  Panthera tigris, and Felis nigripes.
![LYNX](https://s3.us-west-2.amazonaws.com/secure.notion-static.com/1799e45e-d0df-4040-98d5-56171c181265/Untitled.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20220328%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220328T131446Z&X-Amz-Expires=86400&X-Amz-Signature=fc4a219e311265855b8863350987eec4df8e1585c0cfcd69073f753ec62a34df&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22Untitled.png%22&x-id=GetObject)

## Construct genomic windows

### Genome sampling
There are several ways to sample a genome in order to perform a phylogenetic analysis. Just to cite some of examples, you can use SNP-based data, genes (introns, exons, or entire genes), or split the genome into windows of different sizes. For this tutorial, we are going to show you how to split the genome into windows. This method has its problems (like all others), but its major strengths are that it does not rely on genome annotation, and that it allows you to use the whole genome (or its non-repetitive portions) instead of a small fraction of it, such as genes.

Data for this session can be downloaded at:
```
mkdir data
cd data
curl -o LCA.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/LCA.felcat8.chrE1.fasta.gz
curl -o LLY.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/LLY.felcat8.chrE1.fasta.gz
curl -o LRU.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/LRU.felcat8.chrE1.fasta.gz
curl -o PPA.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/PPA.felcat8.chrE1.fasta.gz
curl -o PTI.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/PTI.felcat8.chrE1.fasta.gz
curl -o FNI.felcat8.chrE1.fasta.gz https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-1%20Phylogenomics/data/FNI.felcat8.chrE1.fasta.gz
gunzip *.fasta.gz
sed -i 's/>chrE1:1-61081816/>chrE1/g' FNI.felcat8.chrE1.fasta
```

The initial step will not be shown in this tutorial is to generate the consensus sequence of the genomes. In order for this method to work, the sequence reads from all species must be mapped against the same reference.  To generate the consensus sequence from the mapping files (BAM) you can use [ANGSD](ttp://www.popgen.dk/angsd/index.php/ANGSD) with the "-doFasta 1" parameter. 

### Masking
One important step that is often forgotten is masking low complexity (i.e. repetitive) regions. Those regions usually do not yield reliable alignments and add noise to the results. We are not going to show these steps for this tutorial since they can take quite a while, but it is something you should take into account. Tools like [RepeatMasker](https://www.repeatmasker.org/) can be used to conduct this process.

### Building the genomics windows
The first step is to create an index for the genomes. This is a way that most genomics software use to parse the data and manage computational resources. This can be achieved by running the samtools faidx command.
```
samtools faidx LCA.felcat8.chrE1.fasta
samtools faidx LLY.felcat8.chrE1.fasta
samtools faidx LRU.felcat8.chrE1.fasta
samtools faidx PPA.felcat8.chrE1.fasta
samtools faidx PTI.felcat8.chrE1.fasta
samtools faidx FNI.felcat8.chrE1.fasta
```

If you have more samples and need to speed up the process you can do a loop to create the index of all files with one command.
```
for i in *.fasta; do samtools faidx $i.fasta; done
```

In order to split the genome, you must give the instructions about where to cut the genome. For this part, we are going to use [bedtools](https://bedtools.readthedocs.io/en/latest/), really useful software that manages sequence files (FASTA) and interval files (BED).

The first step is to create a bed file with the genome coordinates that we are interested in. The genome file is a simple text file with the chromosomes or regions and their total size. We can use the index to retrieve this information.
```
cut -f1,2 LCA.felcat8.chrE1.fasta.fai > genome.txt
```

The next step is to run bedtools makewindows, with the genome file, the window size, and the interval between each window, called step size. Window size is total size of the sequences we are going to work with. Step size is the length that we are going to skip between each window. For this tutorial we are using a big step size so we can have fewer windows, but normally use you can apply a much smaller step size or even not use any (in this case, windows would be contiguous).
```
bedtools makewindows -g genome.txt -w 10000 -s 5000000 > 10kb.windows.bed
```
The output looks like this:
```
chrE1	0	10000
chrE1	5000000	5010000
chrE1	10000000	10010000
```
How many lines do we have?

With the bed file that we have created by running the makewindows command, now we can split the genome into windows. This is done by running bedtools getfasta. We use the genome as input, the bed file with the intervals, and set an output.
```
bedtools getfasta -fi LCA.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo LCA.felcat8.chrE1.10kb.fasta
bedtools getfasta -fi LLY.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo LLY.felcat8.chrE1.10kb.fasta
bedtools getfasta -fi LRU.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo LRU.felcat8.chrE1.10kb.fasta
bedtools getfasta -fi PPA.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo PPA.felcat8.chrE1.10kb.fasta
bedtools getfasta -fi PTI.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo PTI.felcat8.chrE1.10kb.fasta
bedtools getfasta -fi FNI.felcat8.chrE1.fasta -bed 10kb.windows.bed -fo FNI.felcat8.chrE1.10kb.fasta
```

To see if everything is ok we can run the following command.
```
grep -c “>” LCA.felctat8.chrE1.fasta
```

The 10kb fasta files have headers formatted this way:
```
>chrE1:0-10000
>chrE1:5000000-5010000
>chrE1:10000000-10010000
```

### Alignment building
The next step is to build the alignments that will be used as inputs for IQTREE (the phylogenetic software – see below). In order to build alignments, we are going to use a custom python script. The logic behind the script is to create files comparing the fasta headers and adding the species name to them.

The script requires that you create a file with the path for the 10kb fasta files. We can achieve this by doing:
```
realpath *.10kb.fasta > dataset.list
```

Now we can run the script. You can download the script to build the alignments using this command:
```
curl -o buildAlignment.py https://raw.githubusercontent.com/henriquevf/buildAlignment/main/buildAlignment.py
python buildAlignment.py dataset.list
```

The script will create a folder with one alignment file for each window containing the five species. Now we are ready to perform the phylogenetic analysis. The fasta header will look like this:
```
>PTI
>PPA
>LCA
>LRU
>LLY
>FNI
```

## Phylogenetic inference based on whole-genome sequences

### Phylogenetic Analysis
Until five to ten years ago, there was a really common discussion about how to use whole-genome sequence data for phylogenetic inference. At that point, the usual way to go was distance-based methods like neighbor-joining. Nowadays, with access to more powerful servers and new algorithms, there is no reason to use distance-based methods (which rely on a simplification of the data) for most projects. The most commonly used software for phylogenetic inference based on molecular data, IQTREE and RAXML employ the maximum likelihood optimality criterion, and can now be run with whole-genome data.

For this tutorial, we are going to use IQTREE. From my experience, IQTREE is more friendly to beginners without losing any of the advanced resources that RAXML has. If you have a small dataset they even have a [web](http://iqtree.cibiv.univie.ac.at/) version.

In order to run the program, we use the following command.
```
iqtree -s chrE1:0-10000.fasta
```
If you need to repeat the analysis add the parameter “-redo”. The program will recognize that you have result files on your folder and will NOT run without this parameter.

The following outputs will be created:
```
chrE1:0-10000.fa.ckp.gz
chrE1:0-10000.fa.log
chrE1:0-10000.fa.model.gz
chrE1:0-10000.fa.bionj
chrE1:0-10000.fa.iqtree
chrE1:0-10000.fa.mldist
chrE1:0-10000.fa.treefile
```
```
Analysis results written to: 
  IQ-TREE report:                chrE1:0-10000.fa.iqtree
  Maximum-likelihood tree:       chrE1:0-10000.fa.treefile
  Likelihood distances:          chrE1:0-10000.fa.mldist
  Screen log file:               chrE1:0-10000.fa.log
```
```
Log-likelihood of the tree: -7857.9437 (s.e. 92.7629)
Unconstrained log-likelihood (without tree): -26429.2435
Number of free parameters (#branches + #model parameters): 11
Akaike information criterion (AIC) score: 15737.8875
Corrected Akaike information criterion (AICc) score: 15737.9139
Bayesian information criterion (BIC) score: 15817.2012

Total tree length (sum of branch lengths): 0.0524
Sum of internal branch lengths: 0.0201 (38.3219% of tree length)

NOTE: Tree is UNROOTED although outgroup taxon 'LCA.felcat8.chrE1.10kb.fasta' is drawn at root

+------------------------------------LCA.felcat8.chrE1.10kb.fasta
|
|     +-----------LLY.felcat8.chrE1.10kb.fasta
+-----|
|     |                                             +-----PPA.felcat8.chrE1.10kb.fasta
|     +---------------------------------------------|
|                                                   +-----PTI.felcat8.chrE1.10kb.fasta
|
+---------------------LRU.felcat8.chrE1.10kb.fasta

Tree in newick format:

(LCA.felcat8.chrE1.10kb.fasta:0.0142199365,(LLY.felcat8.chrE1.10kb.fasta:0.0045896436,(PPA.felcat8.chrE1.10kb.fasta:0.0025620000,PTI.felcat8.chrE1.10kb.fasta:0.0024335088):0.0176021288):0.0024787663,LRU.felcat8.chrE1.10kb.fasta:0.0085145546);

TIME STAMP
----------

Date and time: Fri Mar 11 14:20:27 2022
Total CPU time used: 0.138083 seconds (0h:0m:0s)
Total wall-clock time used: 0.145148 seconds (0h:0m:0s)
```
```
5
LCA.felcat8.chrE1.10kb.fasta 0.0000000 0.0206884 0.0240018 0.0323625 0.0239331 
LLY.felcat8.chrE1.10kb.fasta 0.0206884 0.0000000 0.0135376 0.0244720 0.0274451 
LRU.felcat8.chrE1.10kb.fasta 0.0240018 0.0135376 0.0000000 0.0277762 0.0213774 
PPA.felcat8.chrE1.10kb.fasta 0.0323625 0.0244720 0.0277762 0.0000000 0.0051495 
PTI.felcat8.chrE1.10kb.fasta 0.0239331 0.0274451 0.0213774 0.0051495 0.0000000
```

### Loop
```
for i in *.fa; do iqtree -s $i; done
```

### Bootstrap
A more realistic analysis would include bootstrap to evaluate node support. You can add bootstrap to your trees by running iqtree with the "-bb" parameter. It will do one thousand pseudoreplicates and show the support for each node. 

### Tree parsing
For this step, we will use the software [newick_utils](http://gensoft.pasteur.fr/docs/newick-utils/1.6/nwutils_tutorial.pdf). It is a really powerful tool that can be used directly from the terminal. There are several customization options that are not going to be shown here.
```
nw_display chrE1:0-10000.fa.treefile
```

You can make it look nicer
```
nw_reroot chrE1:0-10000.fasta.treefile FNI | nw_display -w 100 -
```
![tree](https://s3.us-west-2.amazonaws.com/secure.notion-static.com/29eabb64-240f-4af5-be64-4d7b221ebcc3/Untitled.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20220328%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220328T133229Z&X-Amz-Expires=86400&X-Amz-Signature=775b0249e0a5b4907d6ee1ae642aae9bc2c2fc9e5f1c2ca11043d6f9f9b93e63&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22Untitled.png%22&x-id=GetObject)

Another option is to see basic stats about the tree
```
nw_stats chrE1:0-10000.fa.treefile
```
```
Type:	Phylogram
#nodes:	8
#leaves:	5
#dichotomies:	2
#leaf labels:	5
#inner labels:	2
```

This is really handy if you need to compare several trees at once. This can be achieved with the following command. 
```
for i in *.treefile; do nw_stats -f l $i; done
```
```
Phylogram	8	5	2	5	2
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	2
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
Phylogram	8	5	2	5	0
```

A common way to parse data is to count the number of different observed topologies.  This can be done using nw_topology.
```
nw_topology chrE1:0-10000.fa.treefile
```

Now we can stream several commands in order to count how many topologies we have in our dataset.
```
for i in *.treefile; do nw_topology $i; done | nw_order - | uniq -c
```
```
13 (FNI,((LCA,LLY),LRU),(PPA,PTI));
```