# sRNA-spike-in design and analysis
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
To determine the folding structures (i.e. mature.miRNA.seqs.top50percent_folded) and minimum free energies (MFEs) (i.e. mature.miRNA.seqs.top50percent_mfes) of the top 50% expressed mature miRNAs, run the following command:
```shell
python methods.sRNA.spike.in.design.step2.endogenous.miRNAs.py
```
To determine the folding structures (i.e. randomOligoSets/folded/X.fa_folded) and minimum free energies (MFEs) (i.e. randomOligoSets/folded/mfes/X_mfes) of semi-randomly generated sequences from above (i.e. noMatch.fa), run the following command:
```shell
python methods.sRNA.spike.in.design.step2.py
```
To examine distributions of MFEs for each randomly generated set, as well as for the sets specifically used in this manuscript and generated for Supplementary Figure 4 (../supplemental.figure.4/supplemental.figure.4.pdf) run the following: 
```shell
Rscript examine.MFE.distributions.R
```
Please note that the MFE files used to generate the supplemental.figure.4.pdf have been deposited in the randomOligoSets/folded/mfes/ directory, but these will be different each time the above scripts are run because the sequences are randomly generated.


2. methods.data.analysis contains annotations and data analysis scripts to generate data files for small RNA-Seq and mRNA-Seq analyses.
Run the following command in the methods.data.analysis directory with the sample fastq file (e.g. col0_fb1.fa) downloaded from NCBI GEO (Series GSE98553) located within the sample directory (e.g. col0_fb1) to filter raw fastq file. Please note that sample directories and embedded sample.fastq files have to manually created and uploaded, respectively:
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

Now run the following command to convert files generated from the previous step (i.e. allgHits) to individual files containing collections of loci that reads map to (i.e. genomic hit collections) for miRNA gene, tasiRNA gene, transposon and sRNA spike-in loci.
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
To generate files within the data.for.graphs directory with the RPM levels of miRNAs and tasiRNAs coming from each annotated precursor, run the following:
```shell
python getMirAndTasPrecursorVals.py
```
To generate files within each sample directory containing the RPM levels of the sample's sRNA-Seq reads of a specified range as well as the base frequencies of the first read position, run the following:
```shell
python getLengthVals.py col0_leaf1,col0_leaf2,col0_fb1,col0_fb2,d234_fb1,d234_fb2 18 31 all 0
```


3. figure.1 contains files used to generate graphs, as well as example graphs, shown in Figure 1.
Run the following to generate the RPM vs. MPU scatterplots shown in Figure 1 and Supplementary Figure 1:
```shell
Rscript generate.RPM.vs.MPU.scatter.sRNAs.R
```
To generate the density plot shown in Figure 1C, run the following command:
```shell
Rscript generate.miRNA.density.R
```


4. figure.2 contains files used to generate graphs, as well as example graphs, shown in Figure 2.
Run the following to generate violin and scatter plots shown in Figure 2:
```shell
Rscript generate.sRNA.violins.and.miRNA.scatter.R
```


5. figure.3 contains files used to generate graphs, as well as example graphs, shown in Figure 3.
To generate scatter plots shown in Figure 3A and Supplementary Figure 3A-3E, run the following:
```shell
Rscript generate.TPM.vs.MPU.mRNA.R
```
Run the following command to create the one-dimensional strip chart shown in Figure 3B:
```shell
Rscript generate.precursor.strips.R
```
To generate the violin plots showin in Figure 3C, run the following:
```shell
Rscript generate.target.violins.R
```


6. supplemental.figure.1 contains example graphs shown in Supplementary Figure 1.
The example scatterplots were generated in Step 3 above.


7. supplemental.figure.2 contains files used to generate graphs, as well as example graphs, shown in Supplementary Figure 2.
To generate the length distributions shown in Supplementary Figure 2A-2C, run the following command:
```shell
Rscript generate.length.distributions.R
```
Run the following to generate a table (i.e. methods.data.analysis/small.RNA.Seq/data.for.graphs/total.rpm.mlc.table) of total RPM and MPU levels of miRNAs, tasiRNAs and siRNAs in the various samples:
```shell
Rscript generate.total.rpm.mlc.table.R
```
To generate one-dimensional strip charts as shown in Supplementary Figure 2D-2E, run the following command:
```shell
Rscript generate.rpm.and.mpu.strips.R
```


8. supplemental.figure.3 contains example graphs shown in Supplementary Figure 3.
The example scatterplots were generated with generate.TPM.vs.MPU.mRNA.R in Step 5 above.


9. supplemental.figure.4 contains an graph shown in Supplementary Figure 4.
The example violin plot was generated with examine.MFE.distributions.R in Step 1 above.
