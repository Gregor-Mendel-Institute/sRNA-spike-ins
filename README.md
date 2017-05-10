# sRNA-spike-ins
This repository contains scripts used to design and analyze a novel set of small RNA spike-in oligonucleotides (sRNA spike-ins). sRNA spike-ins can be used to normalize high-throughput small RNA sequencing (sRNA-Seq) data across individual experiments, as well as the genome-wide estimation of sRNA:mRNA stoichiometries when used together with mRNA spike-ins. The scripts and files contained within this repository were used to generate all of the data presented in an unpublished manuscript by Lutzmayer, Enugutti and Nodine. The raw data files can be downloaded from the Gene Expression Omnibus (GEO Series GSE98553), and used to generate all data files and graphs according to the procedures outlined below.

Scripts and example graphs are organized into nine main directories:

1. methods.sRNA.spike.in.design contains files used to design and evaluate sRNA spike-ins.

Run the following command in the methods.sRNA.spike.in.design directory to generate 1) a fasta file of the top 50% expressed mature miRNAs (i.e. mature.miRNA.seqs.top50percent.fa) and 2) a fasta file of 21 base sequences with characteristics similar to endogenous miRNAs (i.e. random.fa): 
```shell
python methods.sRNA.spike.in.design.step1.py
```

Use Bowtie to select entries in random.fa generated from above step that do not map to Arabidopsis thaliana (TAIR10) genome, and output to a file called noMatch.fa.

To determine the folding structures (i.e. mature.miRNA.seqs.top50percent_folded) and minimum free energies (MFEs) (i.e. mature.miRNA.seqs.top50percent_mfes) of the top 50% expressed mature miRNAs, run the following command in the methods.sRNA.spike.in.design directory:
```shell
python methods.sRNA.spike.in.design.step2.endogenous.miRNAs.py
```
To determine the folding structures (i.e. randomOligoSets/folded/X.fa_folded) and minimum free energies (MFEs) (i.e. randomOligoSets/folded/mfes/X_mfes) of semi-randomly generated sequences from above (i.e. noMatch.fa), run the following command in the methods.sRNA.spike.in.design directory:
```shell
python methods.sRNA.spike.in.design.step2.py
```
To examine distributions of MFEs for each randomly generated set, as well as for the sets specifically used in this manuscript and generated for Supplementary Figure 4 (../supplemental.figure.4/supplemental.figure.4.pdf) run the examine.MFE.distributions.R script. Please note that the MFE files used to generate the supplemental.figure.4.pdf have been deposited in the randomOligoSets/folded/mfes/ directory, but these will be different each time the above scripts are run.

2. methods.data.analysis contains data analysis pipeline and processed data files for small RNA-Seq and mRNA-Seq analyses.
3. figure.1 contains files used to generate graphs shown in Figure 1.
4. figure.2 contains files used to generate graphs shown in Figure 2.
5. figure.3 contains files used to generate graphs shown in Figure 3. 
6. supplemental.figure.1 contains files used to generate graphs shown in Supplemental Figure 1.
7. supplemental.figure.2 contains files used to generate graphs shown in Supplemental Figure 2.
8. supplemental.figure.3 contains files used to generate graphs shown in Supplemental Figure 3.
9. supplemental.figure.4 contains files used to generate graphs shown in Supplemental Figure 4.
