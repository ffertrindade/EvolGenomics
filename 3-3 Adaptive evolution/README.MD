## 3-3 Adaptive evolution
Today you are going to learn how to perform analyses of natural selection on complete genomes and perform some basic enrichment analyses. We're going to use as example the data from [Li et al., 2019](https://academic.oup.com/mbe/article/36/10/2111/5518928), the same we had used in the [Phylogenomics](https://github.com/ffertrindade/EvolGenomics/tree/main/3-1%20Phylogenomics) class.

Positive selection can be divided into phylogenomic methods, which will be the focus of today’s lesson, and population genomic methods, which require extensive sampling of one or more populations. For phylogenomics-based analysis, we require a whole nuclear genome fasta sequence from three or more species, usually closely related, and the reference genome annotation. If you are working with distantly related species, it is recommended to use *de novo* assemblies with individual annotations. It is a longer pipeline based on orthology tests.

## Investigating adaptive evolution based on genome-wide datasets
For this tutorial, we are going to use the same dataset from day one of module 3. The main difference will be that instead of sampling the genome into windows we are going to retrieve the coding regions (CDS) from chromosome E1 and build a multiple sequence alignment. 

Create a directory for the analysis
```
mkdir selection_class
```
Download the required data 
```
curl -o buildETE3_branchSite.py https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-3%20Adaptive%20evolution/scripts/buildETE3_branchSite.py
curl -o getPositiveGenes.sh https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-3%20Adaptive%20evolution/scripts/getPositiveGenes.sh
curl -o cats_phylogeny.tre https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-3%20Adaptive%20evolution/data/cats_phylogeny.tre
curl -o GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff https://raw.githubusercontent.com/ffertrindade/EvolGenomics/main/3-3%20Adaptive%20evolution/data/GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff
```
The first step is to create indices for the genome fasta files with samtools.
```
samtools faidx LCA.felcat8.chrE1.fasta
samtools faidx LLY.felcat8.chrE1.fasta
samtools faidx LRU.felcat8.chrE1.fasta
samtools faidx PPA.felcat8.chrE1.fasta
samtools faidx PTI.felcat8.chrE1.fasta
```
If you have more samples and need to speed up the process you can do a loop to create the index of all files with one command.
```
for i in *.fasta; do samtools faidx $i.fasta; done
```

Now is the part where it gets different from the phylogenomics tutorial. Instead of creating the windows bed file, we are going to use the annotation file for the domestic cat in gff format. The domestic cat is an important model species used in several projects with wild cats. The annotation used should be from the same genome used as a reference, and also the same version. For this tutorial, we are using an old version of the cat genome and its annotation (felCat8). 

[Gffread](https://github.com/gpertea/gffread) is a small tool that retrieves the coding sequence in various formats. You can add filters to get better results. In our case, it is good since you can remove problematic sequences, that don't have start and stop codons or are out of frame.
```
gffread v0.12.7. Usage:
gffread [-g <genomic_seqs_fasta> | <dir>] [-s <seq_info.fsize>]
[-o <outfile>] [-t <trackname>] [-r [<strand>]<chr>:<start>-<end> [-R]]
[--jmatch <chr>:<start>-<end>] [--no-pseudo]
[-CTVNJMKQAFPGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]
[-j ][--ids <IDs.lst> | --nids <IDs.lst>] [--attrs <attr-list>] [-i <maxintron>]
[--stream] [--bed | --gtf | --tlf] [--table <attrlist>] [--sort-by <ref.lst>]
[<input_gff>]
```
Here's the parameters we are going to use for the tutorial:
```
Filter, convert or cluster GFF/GTF/BED records, extract the sequence of
transcripts (exon or CDS) and more.
By default (i.e. without -O) only transcripts are processed, discarding any
other non-transcript features. Default output is a simplified GFF3 with only
the basic attributes.

Options:
-C   coding only: discard mRNAs that have no CDS features
-V   discard any mRNAs with CDS having in-frame stop codons (requires -g)
-H   for -V option, check and adjust the starting CDS phase
if the original phase leads to a translation with an
in-frame stop codon
-J   discard any mRNAs that either lack initial START codon
or the terminal STOP codon, or have an in-frame stop codon
(i.e. only print mRNAs with a complete CDS)
--no-pseudo: filter out records matching the 'pseudo' keyword
Clustering:
-M/--merge : cluster the input transcripts into loci, discarding
"redundant" transcripts (those with the same exact introns
and fully contained or equal boundaries)
-g   full path to a multi-fasta file with the genomic sequences
for all input mappings, OR a directory with single-fasta files
(one per genomic sequence, with file names matching sequence names)
-x    write a fasta file with spliced CDS for each GFF transcript
-y    write a protein fasta file with the translation of CDS for each record
```
Here is the commands to extract the coding sequences for each species.
```
gffread -g PTI.felcat8.chrE1.fasta -x PTI.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff

gffread -g PPA.felcat8.chrE1.fasta -x PPA.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff

gffread -g LCA.felcat8.chrE1.fasta -x LCA.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff

gffread -g LRU.felcat8.chrE1.fasta -x LRU.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff

gffread -g LLY.felcat8.chrE1.fasta -x LLY.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff

gffread -g FNI.felcat8.chrE1.fasta -x FNI.felcat8.chrE1.cds.fasta -C -V -H -B -J --no-pseudo -M GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff
```

This command will create a fasta file with every CDS that matches our parameters. But differently from the phylogenomics pipeline, the number of CDS will depend on the quality of the genome mapped against the reference. Since we are assuming orthology based on position, the CDS on each species might be different and the filtering done by gffread might remove a different number of genes on each file. How many CDS regions does each species have? 
```
grep -c “>” *.chrE1.cds.fasta

LCA.felcat8.chrE1.cds.fasta:1189
LLY.felcat8.chrE1.cds.fasta:720
LRU.felcat8.chrE1.cds.fasta:1113
PPA.felcat8.chrE1.cds.fasta:1338
PTI.felcat8.chrE1.cds.fasta:1297
FNI.felcat8.chrE1.cds.fasta:1528
```
Now we need to build our alignment. It’s almost the same process as we did for the phylogenomic analyses on day one. In order to make our life easier, we are going to keep only the alignments where we have the five species. What tool would you use for that task? 

First, we create the list file with the CDS
```
realpath *.cds.fasta > cds_alignments.list
buildAlignment.py cds_alignments.list
```
In order to simplify things we need to move some scripts we are going to use further to the windows folders.
```
cp buildETE3.py cds_alignments/windows/
cp cat_phylogeny.tre cds_alignments/windows/
cp getPositiveGenes.sh cds_alignments/windows/
```
How many alignments do we have? 
```
cd cds_alignments/windows/
ls *.fasta | wc -l
```
Now we select only the alignments with the six species. In order to do this, we count how many sequences we on each file. Select the ones that has the six species. Select only the file names and send it to a file. 
```
grep -c ">" *.fasta | grep ":6" | cut -f 1 -d ":" > complete_cds.list
```
We are going to use this list to build the command line to run the selection analyses. We also need the path for the phylogenetic tree for the group. 
```
python buildETE3.py complete_cds.list cat_phylogeny.tre
```

### Running CodeML 
CodeML is part of the [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) software. 

For this tutorial, we are going to use [ETE3](http://etetoolkit.org/). This software is considered a wrapper, a tool that facilitates running another one adding convenient resources. In the case of ETE3 it streamlines the process making it a lot easier to run commands in parallel. It has a command-line interface but also python libraries so you can develop your own scripts around it.

#### Evolutionary models
CodeML is based on a likelihood ratio test, comparing a null hypothesis with another that assumes some kind of selection. The tests are split into three categories: Site-based, Branch-based and Branch-site. 

For this tutorial, we are going to use the site model and branch-site test. Since we are working with a small dataset we do not expect to have a lot of genes under selection. In the real world, we would only do the branch-site test, which is more robust and less prone to false positives. The only problem is that can take more time to run. 

Beyond the alignment we also need a phylogenetic tree. Since we are using the same dataset as the phylogenomics pipeline we could use one the phylogenetic trees we have obtained running IQtree.

#### Branch-site model
The branch-site model is a more robust option that is widely used on positive selection papers. Instead of just comparing the sites between each species, it does that taking into consideration the phylogenetic tree. The pipeline we are going to use here, analyzes the six species in our dataset.  
```
ete3 evol -t cats_phylogeny.tre --alg rna28062.fasta -i rna28062.png -o rna28062 -v3 --models bsA bsA1 --leaves --resume | tee rna28062.stdout.txt
```
As you can see the main difference is the evolutionary model bsA (selection) and bsA1 (relaxed). 

There are several different outputs, but the most important one is printed at the screen. That's the reason we have that final part of the command (tee FILE.stdout.txt). It saves the screen results into a file. The other output files are the alignment file in .png format, which is a great way to visualize the results. And a folder with the original outputs from CodeML.  If your questions requires more twinking of the CodeML parameters, ETE3 has some options to edit the input files. Also, you can edit the output alignment image to better suit your needs. 

That is the output you see on the screen when you run ETE3

We are not covering in this tutorial but you can also run Slr, another selection scan software.
```
Using: /Users/henrique/miniconda3/envs/selection/bin/ete3_apps/bin/Slr
Using: /Users/henrique/miniconda3/envs/selection/bin/ete3_apps/bin/codeml
```
```
Running CodeML/Slr (1 CPUs)
/Users/henrique/miniconda3/envs/selection/lib/python3.6/site-packages/ete3-3.1.2-py3.7.egg/ete3/evol/evoltree.py:540: 
```
```
processing model bsA (bsA~9a0961106af0ebe04e6ab1b4fb043d8a)
marking branches 1
  (FNI #1,((LCA,LLY),LRU),(PPA,PTI));
marking branches 2
  (FNI,((LCA #1,LLY),LRU),(PPA,PTI));
marking branches 3
  (FNI,((LCA,LLY #1),LRU),(PPA,PTI));
marking branches 4
  (FNI,((LCA,LLY),LRU #1),(PPA,PTI));
marking branches 5
  (FNI,((LCA,LLY),LRU),(PPA #1,PTI));
marking branches 6
  (FNI,((LCA,LLY),LRU),(PPA,PTI #1));
```
```
processing model bsA1 (bsA1~9ba11826a26b66cd0bfdf0b567f74acf)
marking branches 1 
  (FNI #1,((LCA,LLY),LRU),(PPA,PTI));
marking branches 2
  (FNI,((LCA #1,LLY),LRU),(PPA,PTI));
marking branches 3
  (FNI,((LCA,LLY #1),LRU),(PPA,PTI));
marking branches 4
  (FNI,((LCA,LLY),LRU #1),(PPA,PTI));
marking branches 5
  (FNI,((LCA,LLY),LRU),(PPA #1,PTI));
marking branches 6
  (FNI,((LCA,LLY),LRU),(PPA,PTI #1));
```
```
LRT
                              Null model |                             Alternative model | p-value
```
```
Null model |                                     Alternative model                               | p-value
   =============================================|===============================================|============
        bsA1.1-1~1f1017a2f939b17f702dba3aae9670d5 |      bsA.1-1~2859df10af171e6b3fe09d562183ab23 | 0.021298*
        bsA1.2-1~07f5abdf541a99cea6a876b2770ad095 |      bsA.2-1~e4ee680f65fe5ce2ad748f18cfd1d71d | 1.000000
        bsA1.3-1~1ce839e8892d3b2a2826e7b106b93268 |      bsA.3-1~5ab11373969e4ff570d430e76bd623df | 1.000000
        bsA1.4-1~e4e837933c4bb721994cc1a31be7da4d |      bsA.4-1~3ef892bd479d0186096fa6f7a178662c | 1.000000
        bsA1.5-1~78645faec0dfea6b6547207de6e057b2 |      bsA.5-1~2b9d87ccfafbc9dbf243c5eb2abca5ee | 1.000000
        bsA1.6-1~4168bca4aeafe68e777f28d1cbc6cf6d |      bsA.6-1~3638b37d8e4159ae57f27b7fec130813 | 1.000000
```
```
SUMMARY BY MODEL
Model bsA.1-1
Marked branches
(FNI #1,((LCA,LLY),LRU),(PPA,PTI));
Model bsA.2-1
Marked branches
(FNI,((LCA #1,LLY),LRU),(PPA,PTI));
Model bsA.3-1
Marked branches
(FNI,((LCA,LLY #1),LRU),(PPA,PTI));
Model bsA.4-1
Marked branches
(FNI,((LCA,LLY),LRU #1),(PPA,PTI));
Model bsA.5-1
Marked branches
(FNI,((LCA,LLY),LRU),(PPA #1,PTI));
Model bsA.6-1
Marked branches
(FNI,((LCA,LLY),LRU),(PPA,PTI #1));
```
```
Marked branches
(FNI #1,((LCA,LLY),LRU),(PPA,PTI));
Sites significantly caracterized
codon position | category
       24   |   Relaxed (probability > 0.99)
Model bsA1.2-1
Marked branches
(FNI,((LCA #1,LLY),LRU),(PPA,PTI));
Sites significantly caracterized
codon position | category
   1-  89   |   Relaxed (probability > 0.99)
Model bsA1.3-1
Marked branches
(FNI,((LCA,LLY #1),LRU),(PPA,PTI));
Sites significantly caracterized
codon position | category
   1-  89   |   Relaxed (probability > 0.99)
Model bsA1.4-1
Marked branches
(FNI,((LCA,LLY),LRU #1),(PPA,PTI));
Sites significantly caracterized
codon position | category
   1-  89   |   Relaxed (probability > 0.99)
Model bsA1.5-1
Marked branches
(FNI,((LCA,LLY),LRU),(PPA #1,PTI));
   1-  89   |   Relaxed (probability > 0.99)
Model bsA1.6-1
Marked branches
(FNI,((LCA,LLY),LRU),(PPA,PTI #1));
Sites significantly caracterized
codon position | category
   1-  89   |   Relaxed (probability > 0.99)
```
We can check the results by running this bash script on the same folder. 
```
bash getPositiveGenes.sh 
```
This script is just a series of bash commands streamlined into a pipeline in order to parse ETE3 ouputs into a table.  If you have a lot of results you can add the species as a filter and have a separate table for each species.
```
grep "bsA.1" .stdout.txt | grep "|" | awk '{print "FNI\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g'
```
And here's our main output. As you can see we don't have a lot of genes under selection. This type of analysis can be run on your machine, but can take a lot of time.  I left this analysis running on my machine for 24 hours, and from about 1600 genes, only about 300 finished. 
```
FNI	rna28508	 0.023561*
FNI	rna28704	 0.021298*
LCA	rna28792	 0.003199**
LCA	rna28793	 0.003199**
PPA	rna28500	 0.015542*
```
Unfortunately, our annotation does not use the gene names. So we need to check the annotation to retrieve the gene names. We pipe the results to a new file with just the sequence names on it. A better way of doing this is creating a master table, so you can merge all this information into one file. 
```
bash getPositiveGenes.sh | cut -f2 > ../../positiveGenes_bs.txt
```
This is what the file looks like
```
rna28508
rna28704
rna28792
rna28793
rna28500
rna28988
rna29308
```
If you want to do manually, you can the command below. First we go back to our root folder. 
```
cd ../../
grep "rna28793" GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff | grep "mRNA”
```
```
chrE1	Gnomon	mRNA	22796935	22807377	.	+	.	ID=rna28793;geneID=gene19902;gene_name=SLFN5
```
To make our life easier in the long run we can add a couple more bash command so we generate a table. This can also be transformed into a bash script.
```
grep -f positiveGenes_bs.txt GCF_000181335.2_Felis_catus_8.0_genomic.chrE1.edited.gff | grep "
mRNA" | cut -f 9 -d$'\t' | cut -f 1,3 -d ";" | sed "s/;/\t/g" | sed "s/ID=//g" | sed "s/gene_name=//g"
```
```
rna28500	LOC101095207
rna28508	LOC101095929
rna28704	LYRM9
rna28792	SLFN5
rna28793	SLFN5
rna28988	STXBP4
rna29308	LOC101087114
```
Here is an example of the alignment from a gene under positive selection. 
![ALIGN](https://s3.us-west-2.amazonaws.com/secure.notion-static.com/a2392a48-141c-4403-9b38-b859f30eb970/Untitled.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20220401%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220401T191008Z&X-Amz-Expires=86400&X-Amz-Signature=eab03ef563a02c74a2ffd545683b013fc67e32b4127d8645dbde48cce5f76613&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22Untitled.png%22&x-id=GetObject)
![ALIGN2](https://s3.us-west-2.amazonaws.com/secure.notion-static.com/3f3efb33-9318-4071-a404-929ad568af60/Untitled.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20220401%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220401T191032Z&X-Amz-Expires=86400&X-Amz-Signature=45e1f8bd5c5ed983957e6fb2a0ed46128bbdfe6b30f6b5b2d3bd7129505eac5a&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22Untitled.png%22&x-id=GetObject)

## Some usefull links and databases and tools for enrichment
During the class, we are going to explore the following tools:
- http://www.webgestalt.org
- https://string-db.org