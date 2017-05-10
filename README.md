# sRNA-spike-ins
This repository contains scripts used to design and analyze a novel set of small RNA spike-in oligonucleotides (sRNA spike-ins). sRNA spike-ins can be used to normalize high-throughput small RNA sequencing (sRNA-Seq) data across individual experiments, as well as the genome-wide estimation of sRNA:mRNA stoichiometries when used together with mRNA spike-ins. The scripts and files contained within this repository were used to generate all of the data presented in an unpublished manuscript by Lutzmayer, Enugutti and Nodine. The raw data files can be downloaded from the NCBI Gene Expression Omnibus (GEO Series GSE98553), and used to generate all data files and graphs according to the procedures outlined below.

Scripts and example graphs are organized into nine main directories in which each set of programs should be run:

1. methods.sRNA.spike.in.design contains files used to design and evaluate sRNA spike-ins.

Run the following command to generate 1) a fasta file of the top 50% expressed mature miRNAs (i.e. mature.miRNA.seqs.top50percent.fa) and 2) a fasta file of 21 base sequences with characteristics similar to endogenous miRNAs (i.e. random.fa): 
```shell
python methods.sRNA.spike.in.design.step1.py
```
Use Bowtie to select entries in random.fa generated from above step that do not map to Arabidopsis thaliana (TAIR10) genome or spike-ins (i.e. genome), and output to a file called noMatch.fa.
```shell
bowtie -f -v 0 --un noMatch.fa genome random.fa match.fa
```
To determine the folding structures (i.e. mature.miRNA.seqs.top50percent_folded) and minimum free energies (MFEs) (i.e. mature.miRNA.seqs.top50percent_mfes) of the top 50% expressed mature miRNAs, run the following command in the methods.sRNA.spike.in.design directory:
```shell
python methods.sRNA.spike.in.design.step2.endogenous.miRNAs.py
```
To determine the folding structures (i.e. randomOligoSets/folded/X.fa_folded) and minimum free energies (MFEs) (i.e. randomOligoSets/folded/mfes/X_mfes) of semi-randomly generated sequences from above (i.e. noMatch.fa), run the following command in the methods.sRNA.spike.in.design directory:
```shell
python methods.sRNA.spike.in.design.step2.py
```
To examine distributions of MFEs for each randomly generated set, as well as for the sets specifically used in this manuscript and generated for Supplementary Figure 4 (../supplemental.figure.4/supplemental.figure.4.pdf) run the following: 
```shell
Rscript examine.MFE.distributions.R
```
Please note that the MFE files used to generate the supplemental.figure.4.pdf have been deposited in the randomOligoSets/folded/mfes/ directory, but these will be different each time the above scripts are run.


2. methods.data.analysis contains annotations and data analysis scripts to generate data files for small RNA-Seq and mRNA-Seq analyses.
Run the following command in the methods.data.analysis directory with sample fastq file (i.e. sample.fastq) downloaded from NCBI GEO (Series GSE98553) located within sample directory to filter raw fastq file. Please note that sample directories and embedded sample.fastq files have to manually created and uploaded, respectively:
```shell
python sRNA_step1_var_qual.py sample 18 75 AGATCGGAAGA no 
```
In the methods.data.analysis directory, run the following command to generate a condensed fasta file where each sequence is represented once and the number of reads per sequence is noted in the entry name:
```shell
python sRNA_step2_v04.py sample 
```
Now use Bowtie to align the fasta file from above (i.e. sample.trimmed.fa) to the Arabidopsis thaliana (TAIR10) genome and ERCC spike-ins (i.e. genome):
```shell
bowtie -f -v 0 -m 100 genome sample/sample.trimmed.fa sample/tags.genome.bwt  
```
Run the following command to convert Bowtie output (i.e. tags.genome.bwt) to tab-delimited file of read attributes including hit-based and reads per million genome-matching read normalizations:
```shell
python bwtToAllHits_v03.py filepath sample total.number.of.genome.matching.reads 
```
Please note that the total.number.of.genome.matching.reads is the number derived from the Bowtie alignment above. Also, the filepath is the one leading up to the directory containing the various sample folders.

Now run the following command to convert file from above step (i.e. allgHits) to individual files containing collections of loci that reads map to (i.e. genomic hit collections) for miRNA gene, tasiRNA gene, transposon and sRNA spike-in loci.
```shell
python make_ghcsAndlDicts_slim.py sample yes yes filepath
```
To generate various files that organize reads that overlap miRNA genes (see comments in script for more details), run the following:
```shell
python miR_finder_v03.py sample filepath
```
To generate various files that organize reads that overlap tasiRNA genes (see comments in script for more details), run the following:
```shell
python tasiR_finder_anno_v03.py sample filepath
```
To generate files within the data.for.graphs directory with the RPM levels of miRNA families, tasiRNA families and 20-22 base or 23-24 base reads contained only within annotated transposons run the following:
```shell
python getMirsSirsAndTas_v02.py
```
To generate files within the data.for.graphs directory with the RPM levels of mature miRNAs, mature tasiRNAs and 20-22 base or 23-24 base reads contained within annotated transposons run the following:
```shell
python getMirsSirsAndTas_v03.py
```
To generate files within the data.for.graphs directory with the levels of miRNAs and tasiRNAs coming from each annotated precursor, run the following:
```shell
python getMirAndTasPrecursorVals.py
```
To generate files within each sample directory containing the levels of miRNAs and tasiRNAs coming from each annotated precursor, run the following:
```shell
python getLengthVals.py col0_leaf1,col0_leaf2,col0_fb1,col0_fb2,d234_fb1,d234_fb2 18 31 all 0
```


3. figure.1 contains files used to generate graphs, as well as example graphs, shown in Figure 1. All of the following commands should be run in the figure.1 directory.
Run the following to generate the RPM vs. MPU scatterplots shown in Figure 1 and Supplementary Figure 1:
```shell
Rscript generate.RPM.vs.MPU.scatter.sRNAs.R
```
To generate the density plot shown in Fig. 1C, run the following command:
```shell
Rscript generate.miRNA.density.R
```


4. figure.2 contains files used to generate graphs shown in Figure 2.
5. figure.3 contains files used to generate graphs shown in Figure 3. 
6. supplemental.figure.1 contains files used to generate graphs shown in Supplemental Figure 1.
7. supplemental.figure.2 contains files used to generate graphs shown in Supplemental Figure 2.
8. supplemental.figure.3 contains files used to generate graphs shown in Supplemental Figure 3.
9. supplemental.figure.4 contains files used to generate graphs shown in Supplemental Figure 4.
