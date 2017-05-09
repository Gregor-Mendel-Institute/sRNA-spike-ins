#generate.miRNA.density.R
#generates density plot for just miRNA fams in col0_fb1 for Fig. 1D

install.packages("extrafont")
library("extrafont")
loadfonts()

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

getVals <- function(sample,type){
  sample_1 = paste(sample,'1',sep='')
  sample_2 = paste(sample,'2',sep='')
  
  data_1 = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample_1,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  data_2 = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample_2,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  
  #now fit linear model and use to predict mpe
  fit_spike_log_1 = lm(log(data_1$mlcs,10) ~ log(data_1$rpm..mean.,10))
  fit_spike_log_2 = lm(log(data_2$mlcs,10) ~ log(data_2$rpm..mean.,10))
  
  #extract individual miR/tasiR families and take average of individual ones with at these 1 rpm in both
  sample_1 = paste(sample,'1','_',type,sep='')
  sample_2 = paste(sample,'2','_',type,sep='')
  vals_1 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample_1,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  vals_2 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample_2,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  
  #select "expressed" small RNAs; i.e. at least 1 rpm
  #also take log10 as this will be used for estimations
  sub_1 = subset(vals_1, rpm >= 1)
  sub_2 = subset(vals_2, rpm >= 1)
  
  sub_1 = log(sub_1,10)
  sub_2 = log(sub_2,10)
  
  sub_1 = as.data.frame(sub_1)
  sub_2 = as.data.frame(sub_2)
  
  sub = merge(sub_1,sub_2, by = "row.names")
  row.names(sub) = sub$Row.names
  names(sub) = c('Row.names','biorep1','biorep2')
  sub = subset(sub, select = c(biorep1,biorep2))
  
  #estimate individual molecules  
  coeffs_1 = coefficients(fit_spike_log_1)
  coeffs_2 = coefficients(fit_spike_log_2)
  
  sub_matrix = cbind(coeffs_1[1] + coeffs_1[2]*sub$biorep1, coeffs_2[1] + coeffs_2[2]*sub$biorep2)   #, mean(c(coeffs_1[1] + coeffs_1[2]*sub$biorep1,coeffs_2[1] + coeffs_2[2]*sub$biorep2)))
  dimnames(sub_matrix)=list(c(row.names(sub)),c("biorep1","biorep2"))
  
  #take exponential values and divid by one million for more meaningful values
  final_matrix = (10 ^ sub_matrix)/1000000
  
  
  return(final_matrix)
}

col0_fb_miR = getVals('col0_fb','miR_fams')
col0_fb_miR_vals = log(rowMeans(col0_fb_miR) * 1000000,10)

pdf("figure.1/fb_miRNA_density_log10.pdf",family='Arial', useDingbats=F)
plot(density(col0_fb_miR_vals), type="h", col="grey", main="miRNA levels in flowers", xlab="Molecules per ug", ylab="Density", las=1)
abline(v=median(col0_fb_miR_vals),lty=2,lwd=2,col='black')
dev.off()