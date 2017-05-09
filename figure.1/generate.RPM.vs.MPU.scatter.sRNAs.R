#generate.RPM.vs.MPU.scatter.sRNAs.R
#creates scatter-plots of small RNA spike-in RPM vs. MPU levels as shown in Fig. 1B and Supplementary Fig. 1A-1E.

install.packages("extrafont")
library("extrafont")
#font_import()
loadfonts()

sample = 'col0_fb1'     #sample can be 'col0_fb1', 'col0_fb2', 'col0_leaf1', 'col0_leaf2', 'd234_fb1' or 'd234_fb2'
denominator = 'ug'

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)

#log10-transform, determine correlation and fit linear model on individual datasets
data_log10 = log(data,10)  
cor_log10 = cor(data_log10$rpm..mean., data_log10$mlcs, method=c("pearson"))
cor_log10
cor.test(data_log10$rpm..mean., data_log10$mlcs, method=c("pearson"))
fit_log10 = lm(data_log10$mlcs ~ data_log10$rpm..mean.)

cexMag = 2.5

#now plot individual graphs with log10  
##col0_fb1 goes to figure.1/ and rest go to supplemental.figure.1/ directory
plotFile = paste("figure.1/rpm_vs_mol_",sample,"_log10.pdf", sep = "")
pdf(plotFile, family="Arial", useDingbats=FALSE)
plot(data_log10$rpm..mean., data_log10$mlcs, pch=20, cex=2, cex.lab=cexMag, cex.axis=cexMag, cex.main=cexMag * 0.8571429, las=1, xlab=expression("Reads per million (log"[10]*")"), ylab=expression("Molecules per ug (log"[10]*")"))#mgp=c(3,1.5,0))
text(min(data_log10$rpm..mean.), max(data_log10$mlcs) - 0.2, sample,pos=4, cex=cexMag, font=2)#, adj=10)
text(min(data_log10$rpm..mean.), max(data_log10$mlcs) - 0.5, paste("R =",round(cor_log10, digits = 3)),pos=4, cex=cexMag)#, adj=10)
abline(fit_log10, col="black", lty=2, lwd=1)
dev.off()
 
