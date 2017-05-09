#generate.TPM.vs.MPU.mRNA.R
#creates scatter-plots of mRNA ERCC spike-in RPM vs. MPU levels as shown in Fig. 3A and Supplementary Fig. 3A-3E.

#install.packages("extrafont")
library("extrafont")
#font_import()
loadfonts()

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages('vioplot')
library("vioplot")

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

master = read.table("methods.data.analysis/master_genes_Araport.txt", header=TRUE, sep='\t', row.names=1, strip.white=T)

Na = 6.02214179 * 10e23 * 1e-18

plot.ERCCs <- function(sample,mix,outFile){
  scaling_factor = 0.0002 #scaling_factor was the same for all samples used in this study; scaling_factor = (ul of ERCC mix added) * (1/dilution of total ug that spike-in was added to) * (1/ERCC dilution) * (1/total ug)
                                                                                                          # = (1) * (1/50) * (1/200) * (1/0.5)

  #Analysis.txt was downloaded from LifeTech's website in July, 2013
  spike = read.table("methods.data.analysis/Analysis.txt", header=TRUE, sep="\t", row.names=2, strip.white=T)
  spike = subset(spike, select=paste("concentration.in.Mix.",mix,"..attomoles.ul.",sep=""))
  colnames(spike) = "molecules"
  
  #multiply by Avagadro's number and scaling factor to determine number of molecules for each spike-in added
  spike = spike * Na * scaling_factor
  
  #select ERCCs with values >= 1 TPM
  ERCC_all = merge(master,spike, by = "row.names")
  row.names(ERCC_all) = ERCC_all$Row.names
  ERCC_exp = subset(ERCC_all, ERCC_all[,sample] >= 1, select = c(sample,'molecules'))
  colnames(ERCC_exp) = c("TPM","molecules")
  
  fit_log = lm(log(ERCC_exp$molecules,10) ~ log(ERCC_exp$TPM,10))
  
  #plot graph
  pdf(outFile, family="Arial", useDingbats=FALSE)
  plot(log(ERCC_exp$TPM,10), log(ERCC_exp$molecules,10), pch=20, cex=2, las=1, main=sample, xlab=c("Transcripts per million (log10)"), ylab=c("Molecules per ug (log10)"))
  text(0, max(log(ERCC_exp$molecules,10)) - 0.25, paste("Pearson's R =",round(cor(log(ERCC_exp$TPM,10),log(ERCC_exp$molecules,10))[1],digits=3)),pos=4)#, adj=10)
  abline(fit_log, col="black", lty=2, lwd=2)
  dev.off()
}

plot.ERCCs('col0_fb1',1,'figure.3/col0_fb1.ERCC.TPM.vs.MPU.scatter.pdf')
plot.ERCCs('col0_fb2',1,'supplemental.figure.3/col0_fb2.ERCC.TPM.vs.MPU.scatter.pdf')
plot.ERCCs('col0_leaf1',1,'supplemental.figure.3/col0_leaf1.ERCC.TPM.vs.MPU.scatter.pdf')
plot.ERCCs('col0_leaf2',1,'supplemental.figure.3/col0_leaf2.ERCC.TPM.vs.MPU.scatter.pdf')
plot.ERCCs('d234_fb1',2,'supplemental.figure.3/d234_fb1.ERCC.TPM.vs.MPU.scatter.pdf')
plot.ERCCs('d234_fb2',2,'supplemental.figure.3/d234_fb2.ERCC.TPM.vs.MPU.scatter.pdf')


