# sRNA-spike-ins
This repository contains scripts used to design and analyze a novel set of small RNA spike-in oligonucleotides (sRNA spike-ins). sRNA spike-ins can be used to normalize high-throughput small RNA sequencing (sRNA-Seq) data across individual experiments, as well as the genome-wide estimation of sRNA:mRNA stoichiometries when used together with mRNA spike-ins. The scripts and files contained within this repository were used to generate all of the data presented in an unpublished manuscript by Lutzmayer, Enugutti and Nodine. The raw data files can be downloaded from the Gene Expression Omnibus (GEO Series GSE98553), and used to generate all data files and graphs according to the procedures outlined below.

Scripts and example graphs are organized into nine main directories:

1. methods.sRNA.spike.in.design contains files used to design and evaluate sRNA spike-ins.

To generate a random set of 21 base sequences with characteristics similar to endogenous miRNAs, run this command in the methods.sRNA.spike.in.design directory: 
```shell
python methods.sRNA.spike.in.design.step1.py
```
2. methods.data.analysis contains data analysis pipeline and processed data files for small RNA-Seq and mRNA-Seq analyses.
3. figure.1 contains files used to generate graphs shown in Figure 1.
4. figure.2 contains files used to generate graphs shown in Figure 2.
5. figure.3 contains files used to generate graphs shown in Figure 3. 
6. supplemental.figure.1 contains files used to generate graphs shown in Supplemental Figure 1.
7. supplemental.figure.2 contains files used to generate graphs shown in Supplemental Figure 2.
8. supplemental.figure.3 contains files used to generate graphs shown in Supplemental Figure 3.
9. supplemental.figure.4 contains files used to generate graphs shown in Supplemental Figure 4.
