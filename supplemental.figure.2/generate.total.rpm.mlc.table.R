#generate.total.rpm.mlc.table.R
#calculates total levels of miRNAs, tasiRNAs, 20-22 nt siRNAs and 23-24 nt siRNAs in RPM units, converts these to MPU units and prints total.number of various sRNA classes and prints to "total.rpm.mlc.table" file in generate.rpm.and.mpu.strips.R program used for Supplementary Fig. 2D-2E

#install.packages("extrafont")
library("extrafont")
#font_import()
loadfonts()

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

getVals <- function(sample){
  #select data from premade dose-response tables
  data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)

  siR_24 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample,"_siR_24_v2",sep=""), header=TRUE, sep="\t", row.names=NULL, strip.white=T)
  siR_21 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample,"_siR_21_v2",sep=""), header=TRUE, sep="\t", row.names=NULL, strip.white=T)
  miR = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample,"_miR_mature",sep=""), header=TRUE, sep="\t", row.names=NULL, strip.white=T)
  tasiR = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample,"_tasiR_mature",sep=""), header=TRUE, sep="\t", row.names=NULL, strip.white=T)
  
  #now fit linear model and use to predict mpe
  fit_log_sRNA = lm(log(data$mlcs,10) ~ log(data$rpm..mean.,10))

  coeffs = coefficients(fit_log_sRNA)

  #sub_matrix = cbind(sum(miR$rpm), coeffs[1] + coeffs[2]*sum(miR$rpm), sum(siR_24$rpm), coeffs[1] + coeffs[2]*sum(siR_24$rpm), sum(siR_21$rpm), coeffs[1] + coeffs[2]*sum(siR_21$rpm), sum(tasiR$rpm), coeffs[1] + coeffs[2]*sum(tasiR$rpm))
  sub_matrix = cbind(sum(miR$rpm), 10^(coeffs[1] + coeffs[2]*log(sum(miR$rpm),10)),
                     sum(siR_24$rpm), 10^(coeffs[1] + coeffs[2]*log(sum(siR_24$rpm),10)),
                     sum(siR_21$rpm), 10^(coeffs[1] + coeffs[2]*log(sum(siR_21$rpm),10)),
                     sum(tasiR$rpm), 10^(coeffs[1] + coeffs[2]*log(sum(tasiR$rpm),10)))
                     
  dimnames(sub_matrix)=list(c(sample),c("miR_rpm","miR_mlc","siR_24_rpm","siR_24_mlc","siR_21_rpm","siR_21_mlc","tasiR_rpm","tasiR_mlc"))

  return(sub_matrix)
}

fb1_sub = getVals("col0_fb1")
fb2_sub = getVals("col0_fb2")
leaf1_sub = getVals("col0_leaf1")
leaf2_sub = getVals("col0_leaf2")
d234_fb1_sub = getVals("d234_fb1")
d234_fb2_sub = getVals("d234_fb2")

total = rbind(fb1_sub, fb2_sub, leaf1_sub, leaf2_sub, d234_fb1_sub, d234_fb2_sub)

write.table(total,"methods.data.analysis/small.RNA.Seq/data.for.graphs/total.rpm.mlc.table", sep = '\t', quote=F)
