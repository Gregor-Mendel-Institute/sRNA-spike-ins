#generate.RPM.vs.MPU.scatter.sRNAs.R
#generates violin and scatter plots shown in Fig. 2

#install.packages("extrafont")
library("extrafont")
loadfonts()

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages('vioplot')
library("vioplot")

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

getVals_mpu <- function(sample,type){
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

  sub_matrix = data.frame(rowMeans(cbind(coeffs_1[1] + coeffs_1[2]*sub$biorep1, coeffs_2[1] + coeffs_2[2]*sub$biorep2)))   
  dimnames(sub_matrix)=list(c(row.names(sub)),c("MPU"))

  #take exponential values to generate expected MPU
  final_matrix = (10 ^ sub_matrix)
  
  return(final_matrix)
}

getVals_rpm <- function(sample,type){
  #extract individual miR/tasiR families and take average of individual ones with at these 1 rpm in both
  sample_1 = paste(sample,'1','_',type,sep='')
  sample_2 = paste(sample,'2','_',type,sep='')
  vals_1 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample_1,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  vals_2 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample_2,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  
  #select "expressed" small RNAs; i.e. at least 1 rpm
  sub_1 = subset(vals_1, rpm >= 1)
  sub_2 = subset(vals_2, rpm >= 1)

  sub_1 = as.data.frame(sub_1)
  sub_2 = as.data.frame(sub_2)
  
  sub = merge(sub_1,sub_2, by = "row.names")
  row.names(sub) = sub$Row.names
  names(sub) = c('Row.names','biorep1','biorep2')
  sub = subset(sub, select = c(biorep1,biorep2))

  sub_matrix = data.frame(rowMeans(sub))   
  dimnames(sub_matrix)=list(c(row.names(sub)),c("RPM"))
    
  return(sub_matrix)
}

getVals_sig <- function(tissue1,tissue2){
  sample1_1 = paste(tissue1,'1',sep='')
  sample1_2 = paste(tissue1,'2',sep='')

  sample2_1 = paste(tissue2,'1',sep='')
  sample2_2 = paste(tissue2,'2',sep='')
  
  sample1_1_data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample1_1,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  sample1_2_data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample1_2,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  sample2_1_data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample2_1,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  sample2_2_data = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample2_2,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)
  
  #now fit linear model and use to predict mpe
  fit_spike_log1_1 = lm(log(sample1_1_data$mlcs,10) ~ log(sample1_1_data$rpm..mean.,10))
  fit_spike_log1_2 = lm(log(sample1_2_data$mlcs,10) ~ log(sample1_2_data$rpm..mean.,10))
  fit_spike_log2_1 = lm(log(sample2_1_data$mlcs,10) ~ log(sample2_1_data$rpm..mean.,10))
  fit_spike_log2_2 = lm(log(sample2_2_data$mlcs,10) ~ log(sample2_2_data$rpm..mean.,10))
    
  #extract individual miR families and take average of individual ones with at these 1 rpm in both
  sample1_1 = paste(sample1_1,'_miR_fams',sep='')
  sample1_2 = paste(sample1_2,'_miR_fams',sep='')
  sample2_1 = paste(sample2_1,'_miR_fams',sep='')
  sample2_2 = paste(sample2_2,'_miR_fams',sep='')
  
  vals1_1 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample1_1,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  vals1_2 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample1_2,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  vals2_1 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample2_1,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  vals2_2 = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample2_2,sep=""), header=TRUE, row.names=1, sep="\t", strip.white=T, quote="")
  
  #select "expressed" small RNAs; i.e. at least 1 rpm
  sub1_1 = subset(vals1_1, rpm >= 1)
  sub1_2 = subset(vals1_2, rpm >= 1)
  sub2_1 = subset(vals2_1, rpm >= 1)
  sub2_2 = subset(vals2_2, rpm >= 1)
  
  sub1_1 = as.data.frame(sub1_1)
  sub1_2 = as.data.frame(sub1_2)
  sub2_1 = as.data.frame(sub2_1)
  sub2_2 = as.data.frame(sub2_2)
  
  sub = merge(sub1_1,sub1_2, by = "row.names")
  row.names(sub) = sub$Row.names
  sub = subset(sub, select = c(2:3))
  names(sub) = c(tissue1_1,tissue1_2)

  sub = merge(sub,sub2_1, by = "row.names")
  row.names(sub) = sub$Row.names
  sub = subset(sub, select = c(2:4))
  names(sub) = c(tissue1_1,tissue1_2,tissue2_1)
  
  sub = merge(sub,sub2_2, by = "row.names")
  row.names(sub) = sub$Row.names
  sub = subset(sub, select = c(2:5))
  names(sub) = c(tissue1_1,tissue1_2,tissue2_1,tissue2_2)  
  
  #estimate individual molecules  
  coeffs1_1 = coefficients(fit_spike_log1_1)
  coeffs1_2 = coefficients(fit_spike_log1_2)
  coeffs2_1 = coefficients(fit_spike_log2_1)
  coeffs2_2 = coefficients(fit_spike_log2_2)

  all = data.frame(cbind(sub[,1],sub[,2],sub[,3],sub[,4],10^(coeffs1_1[1] + coeffs1_1[2]*log(sub[,1],10)),10^(coeffs1_2[1] + coeffs1_2[2]*log(sub[,2],10)),10^(coeffs2_1[1] + coeffs2_1[2]*log(sub[,3],10)),10^(coeffs2_2[1] + coeffs2_2[2]*log(sub[,4],10))))
  dimnames(all)=list(c(row.names(sub)),c("tissue1_rep1_RPM","tissue1_rep2_RPM","tissue2_rep1_RPM","tissue2_rep2_RPM","tissue1_rep1_MPU","tissue1_rep2_MPU","tissue2_rep1_MPU","tissue2_rep2_MPU"))
  
  #next select miRNAs with significant (t.test) and 2-fold differences for plotting and counting
  ##rpm
  all_sig_rpm = apply(all, 1, function(x) t.test(c(x[1],x[2]),c(x[3],x[4]))$p.value)

  #combine with all from above
  all_sig_rpm_ext = merge(all,all_sig_rpm, by = "row.names")
  row.names(all_sig_rpm_ext) = all_sig_rpm_ext$Row.names
  all_sig_rpm_plus = subset(all_sig_rpm_ext, ((tissue2_rep1_RPM + tissue2_rep2_RPM)/2)/((tissue1_rep1_RPM + tissue1_rep2_RPM)/2) >= 2 & y <= 0.05, select = c(2:5,10))
  all_sig_rpm_minus = subset(all_sig_rpm_ext, ((tissue2_rep1_RPM + tissue2_rep2_RPM)/2)/((tissue1_rep1_RPM + tissue1_rep2_RPM)/2) <= 0.5 & y <= 0.05, select = c(2:5,10))
  
  #take mean of each of two bioreps
  #increased transcripts in tissue 2 vs tissue 1
  sig_rpm_plus = data.frame(t(apply(all_sig_rpm_plus, 1, function(x) matrix(c(mean(c(x[1],x[2])),mean(c(x[3],x[4])), x[5]), byrow=TRUE))))

  if (length(sig_rpm_plus) > 1) {
  rownames(sig_rpm_plus) = rownames(all_sig_rpm_plus)
  names(sig_rpm_plus) = c(tissue1,tissue2,'P')
  }
  
  #decreased transcripts in tissue 2 vs tissue 1
  sig_rpm_minus = data.frame(t(apply(all_sig_rpm_minus, 1, function(x) matrix(c(mean(c(x[1],x[2])),mean(c(x[3],x[4])), x[5]), byrow=TRUE))))
  
  if (length(sig_rpm_minus) > 1) {
  rownames(sig_rpm_minus) = rownames(all_sig_rpm_minus)
  names(sig_rpm_minus) = c(tissue1,tissue2,'P')
  }
  
  ##mpu
  all_sig_mpu = apply(all, 1, function(x) t.test(c(x[5],x[6]),c(x[7],x[8]))$p.value)
  
  #combine with all from above
  all_sig_mpu_ext = merge(all,all_sig_mpu, by = "row.names")
  row.names(all_sig_mpu_ext) = all_sig_mpu_ext$Row.names
  
  all_sig_mpu_plus = subset(all_sig_mpu_ext, ((tissue2_rep1_MPU + tissue2_rep2_MPU)/2)/((tissue1_rep1_MPU + tissue1_rep2_MPU)/2) >= 2 & y < 0.05, select = c(6:10))
  all_sig_mpu_minus = subset(all_sig_mpu_ext, ((tissue2_rep1_MPU + tissue2_rep2_MPU)/2)/((tissue1_rep1_MPU + tissue1_rep2_MPU)/2) <= 0.5 & y < 0.05, select = c(6:10))
  
  #take mean of each of two bioreps
  #increased transcripts in tissue 2 vs tissue 1
  sig_mpu_plus = data.frame(t(apply(all_sig_mpu_plus, 1, function(x) matrix(c(mean(c(log(x[1],10),log(x[2],10))),mean(c(log(x[3],10),log(x[4],10))), x[5]), byrow=TRUE))))
  
  if (length(sig_mpu_plus) > 1) {
  row.names(sig_mpu_plus) = rownames(all_sig_mpu_plus)
  names(sig_mpu_plus) = c(tissue1,tissue2,'P')
  }
  #decreased transcripts in tissue 2 vs tissue 1
  sig_mpu_minus = data.frame(t(apply(all_sig_mpu_minus, 1, function(x) matrix(c(mean(c(log(x[1],10),log(x[2],10))),mean(c(log(x[3],10),log(x[4],10))), x[5]), byrow=TRUE))))
  
  if (length(sig_mpu_minus) > 1) {
  rownames(sig_mpu_minus) = rownames(all_sig_mpu_minus)
  names(sig_mpu_minus) = c(tissue1,tissue2,'P')
  }
  
  return(list(sig_rpm_plus,sig_rpm_minus,sig_mpu_plus,sig_mpu_minus))
}


#extract subsets for individual plotting
##MPU
col0_fb_miR_mpu = getVals_mpu('col0_fb','miR_fams')
col0_leaf_miR_mpu = getVals_mpu('col0_leaf','miR_fams')
d234_fb_miR_mpu = getVals_mpu('d234_fb','miR_fams')

col0_fb_tasiR_mpu = getVals_mpu('col0_fb','tasiR_fams')
col0_leaf_tasiR_mpu = getVals_mpu('col0_leaf','tasiR_fams')
d234_fb_tasiR_mpu = getVals_mpu('d234_fb','tasiR_fams')

col0_fb_siR_21_mpu = getVals_mpu('col0_fb','siR_21')
col0_leaf_siR_21_mpu = getVals_mpu('col0_leaf','siR_21')
d234_fb_siR_21_mpu = getVals_mpu('d234_fb','siR_21')

col0_fb_siR_24_mpu = getVals_mpu('col0_fb','siR_24')
col0_leaf_siR_24_mpu = getVals_mpu('col0_leaf','siR_24')
d234_fb_siR_24_mpu = getVals_mpu('d234_fb','siR_24')

##RPM
col0_fb_miR_rpm = getVals_rpm('col0_fb','miR_fams')
col0_leaf_miR_rpm = getVals_rpm('col0_leaf','miR_fams')
d234_fb_miR_rpm = getVals_rpm('d234_fb','miR_fams')

col0_fb_tasiR_rpm = getVals_rpm('col0_fb','tasiR_fams')
col0_leaf_tasiR_rpm = getVals_rpm('col0_leaf','tasiR_fams')
d234_fb_tasiR_rpm = getVals_rpm('d234_fb','tasiR_fams')

col0_fb_siR_21_rpm = getVals_rpm('col0_fb','siR_21')
col0_leaf_siR_21_rpm = getVals_rpm('col0_leaf','siR_21')
d234_fb_siR_21_rpm = getVals_rpm('d234_fb','siR_21')

col0_fb_siR_24_rpm = getVals_rpm('col0_fb','siR_24')
col0_leaf_siR_24_rpm = getVals_rpm('col0_leaf','siR_24')
d234_fb_siR_24_rpm = getVals_rpm('d234_fb','siR_24')


#merge by common row names
##MPU
###miRNA
miR_merge_mpu = merge(col0_leaf_miR_mpu,col0_fb_miR_mpu, by = "row.names")
row.names(miR_merge_mpu) = miR_merge_mpu$Row.names
miR_merge_mpu = subset(miR_merge_mpu, select = c(2:3))
names(miR_merge_mpu) = c('MPU_leaf','MPU_fb')

miR_merge_mpu = merge(miR_merge_mpu,d234_fb_miR_mpu, by = "row.names")
row.names(miR_merge_mpu) = miR_merge_mpu$Row.names
miR_merge_mpu = subset(miR_merge_mpu, select = c(2:4))
names(miR_merge_mpu) = c('MPU_leaf','MPU_fb','MPU_d234')

###tasiRNA
tasiR_merge_mpu = merge(col0_leaf_tasiR_mpu,col0_fb_tasiR_mpu, by = "row.names")
row.names(tasiR_merge_mpu) = tasiR_merge_mpu$Row.names
tasiR_merge_mpu = subset(tasiR_merge_mpu, select = c(2:3))
names(tasiR_merge_mpu) = c('MPU_leaf','MPU_fb')

tasiR_merge_mpu = merge(tasiR_merge_mpu,d234_fb_tasiR_mpu, by = "row.names")
row.names(tasiR_merge_mpu) = tasiR_merge_mpu$Row.names
tasiR_merge_mpu = subset(tasiR_merge_mpu, select = c(2:4))
names(tasiR_merge_mpu) = c('MPU_leaf','MPU_fb','MPU_d234')

###siR_21
siR_21_merge_mpu = merge(col0_leaf_siR_21_mpu,col0_fb_siR_21_mpu, by = "row.names")
row.names(siR_21_merge_mpu) = siR_21_merge_mpu$Row.names
siR_21_merge_mpu = subset(siR_21_merge_mpu, select = c(2:3))
names(siR_21_merge_mpu) = c('MPU_leaf','MPU_fb')

siR_21_merge_mpu = merge(siR_21_merge_mpu,d234_fb_siR_21_mpu, by = "row.names")
row.names(siR_21_merge_mpu) = siR_21_merge_mpu$Row.names
siR_21_merge_mpu = subset(siR_21_merge_mpu, select = c(2:4))
names(siR_21_merge_mpu) = c('MPU_leaf','MPU_fb','MPU_d234')

###siR_24
siR_24_merge_mpu = merge(col0_leaf_siR_24_mpu,col0_fb_siR_24_mpu, by = "row.names")
row.names(siR_24_merge_mpu) = siR_24_merge_mpu$Row.names
siR_24_merge_mpu = subset(siR_24_merge_mpu, select = c(2:3))
names(siR_24_merge_mpu) = c('MPU_leaf','MPU_fb')

siR_24_merge_mpu = merge(siR_24_merge_mpu,d234_fb_siR_24_mpu, by = "row.names")
row.names(siR_24_merge_mpu) = siR_24_merge_mpu$Row.names
siR_24_merge_mpu = subset(siR_24_merge_mpu, select = c(2:4))
names(siR_24_merge_mpu) = c('MPU_leaf','MPU_fb','MPU_d234')


##RPM
###miRNA
miR_merge_rpm = merge(col0_leaf_miR_rpm,col0_fb_miR_rpm, by = "row.names")
row.names(miR_merge_rpm) = miR_merge_rpm$Row.names
miR_merge_rpm = subset(miR_merge_rpm, select = c(2:3))
names(miR_merge_rpm) = c('RPM_leaf','RPM_fb')

miR_merge_rpm = merge(miR_merge_rpm,d234_fb_miR_rpm, by = "row.names")
row.names(miR_merge_rpm) = miR_merge_rpm$Row.names
miR_merge_rpm = subset(miR_merge_rpm, select = c(2:4))
names(miR_merge_rpm) = c('RPM_leaf','RPM_fb','RPM_d234')

###tasiRNA
tasiR_merge_rpm = merge(col0_leaf_tasiR_rpm,col0_fb_tasiR_rpm, by = "row.names")
row.names(tasiR_merge_rpm) = tasiR_merge_rpm$Row.names
tasiR_merge_rpm = subset(tasiR_merge_rpm, select = c(2:3))
names(tasiR_merge_rpm) = c('RPM_leaf','RPM_fb')

tasiR_merge_rpm = merge(tasiR_merge_rpm,d234_fb_tasiR_rpm, by = "row.names")
row.names(tasiR_merge_rpm) = tasiR_merge_rpm$Row.names
tasiR_merge_rpm = subset(tasiR_merge_rpm, select = c(2:4))
names(tasiR_merge_rpm) = c('RPM_leaf','RPM_fb','RPM_d234')

###siR_21
siR_21_merge_rpm = merge(col0_leaf_siR_21_rpm,col0_fb_siR_21_rpm, by = "row.names")
row.names(siR_21_merge_rpm) = siR_21_merge_rpm$Row.names
siR_21_merge_rpm = subset(siR_21_merge_rpm, select = c(2:3))
names(siR_21_merge_rpm) = c('RPM_leaf','RPM_fb')

siR_21_merge_rpm = merge(siR_21_merge_rpm,d234_fb_siR_21_rpm, by = "row.names")
row.names(siR_21_merge_rpm) = siR_21_merge_rpm$Row.names
siR_21_merge_rpm = subset(siR_21_merge_rpm, select = c(2:4))
names(siR_21_merge_rpm) = c('RPM_leaf','RPM_fb','RPM_d234')

###siR_24
siR_24_merge_rpm = merge(col0_leaf_siR_24_rpm,col0_fb_siR_24_rpm, by = "row.names")
row.names(siR_24_merge_rpm) = siR_24_merge_rpm$Row.names
siR_24_merge_rpm = subset(siR_24_merge_rpm, select = c(2:3))
names(siR_24_merge_rpm) = c('RPM_leaf','RPM_fb')

siR_24_merge_rpm = merge(siR_24_merge_rpm,d234_fb_siR_24_rpm, by = "row.names")
row.names(siR_24_merge_rpm) = siR_24_merge_rpm$Row.names
siR_24_merge_rpm = subset(siR_24_merge_rpm, select = c(2:4))
names(siR_24_merge_rpm) = c('RPM_leaf','RPM_fb','RPM_d234')

#plot rpm-based violin plots
##convert to dataframes to log10 for plotting
miR_merge_rpm_log = log(miR_merge_rpm,10)
tasiR_merge_rpm_log = log(tasiR_merge_rpm,10)
siR_21_merge_rpm_log = log(siR_21_merge_rpm,10)
siR_24_merge_rpm_log = log(siR_24_merge_rpm,10)

miR_merge_mpu_log = log(miR_merge_mpu,10)
tasiR_merge_mpu_log = log(tasiR_merge_mpu,10)
siR_21_merge_mpu_log = log(siR_21_merge_mpu,10)
siR_24_merge_mpu_log = log(siR_24_merge_mpu,10)


name_v = c("Col-0\nleaves","Col-0\nflowers","dcl234\nflowers")
col_v <- brewer.pal(4,"Paired")
col_v = c(col_v[4],col_v[3],col_v[1])

pdf("figure.2/sRNA_violin_rpm.pdf", family='Arial', useDingbats=F)
par(mfrow=c(2,2), las=1)
vioplot(miR_merge_rpm_log$RPM_leaf,miR_merge_rpm_log$RPM_fb,miR_merge_rpm_log$RPM_d234, ylim=c(0,4.7), col='white', names=name_v)
vioplot(miR_merge_rpm_log$RPM_leaf,col=col_v[1], at=1, add=T)
vioplot(miR_merge_rpm_log$RPM_fb,col=col_v[2], at=2, add=T)
vioplot(miR_merge_rpm_log$RPM_d234,col=col_v[3], at=3, add=T)
text(0.4,4.6,labels=c("miRNAs"), pos=4)
abline(h=median(miR_merge_rpm_log$RPM_fb), lty=2, lwd=1, col="grey")
#tasiRNAs
vioplot(tasiR_merge_rpm_log$RPM_leaf,tasiR_merge_rpm_log$RPM_fb,tasiR_merge_rpm_log$RPM_d234, ylim=c(0,4.7), col='white', names=name_v)
vioplot(tasiR_merge_rpm_log$RPM_leaf,col=col_v[1], at=1, add=T)
vioplot(tasiR_merge_rpm_log$RPM_fb,col=col_v[2], at=2, add=T)
vioplot(tasiR_merge_rpm_log$RPM_d234,col=col_v[3], at=3, add=T)
text(0.4,4.6,labels=c("tasiRNAs"), pos=4)
abline(h=median(tasiR_merge_rpm_log$RPM_fb), lty=2, lwd=1, col="grey")
#20-22 nt siRNAs
vioplot(siR_21_merge_rpm_log$RPM_leaf,siR_21_merge_rpm_log$RPM_fb,siR_21_merge_rpm_log$RPM_d234, ylim=c(0,4.7), col='white', names=name_v)
vioplot(siR_21_merge_rpm_log$RPM_leaf,col=col_v[1], at=1, add=T)
vioplot(siR_21_merge_rpm_log$RPM_fb,col=col_v[2], at=2, add=T)
vioplot(siR_21_merge_rpm_log$RPM_d234,col=col_v[3], at=3, add=T)
text(0.4,4.6,labels=c("20-22 nt siRNAs"), pos=4)
abline(h=median(siR_21_merge_rpm_log$RPM_fb), lty=2, lwd=1, col="grey")
#23-24 nt siRNAs
vioplot(siR_24_merge_rpm_log$RPM_leaf,siR_24_merge_rpm_log$RPM_fb,siR_24_merge_rpm_log$RPM_d234, ylim=c(0,4.7), col='white', names=name_v)
vioplot(siR_24_merge_rpm_log$RPM_leaf,col=col_v[1], at=1, add=T)
vioplot(siR_24_merge_rpm_log$RPM_fb,col=col_v[2], at=2, add=T)
vioplot(siR_24_merge_rpm_log$RPM_d234,col=col_v[3], at=3, add=T)
text(0.4,4.6,labels=c("23-24 nt siRNAs"), pos=4)
abline(h=median(siR_24_merge_rpm_log$RPM_fb), lty=2, lwd=1, col="grey")
dev.off()

ks.test(miR_merge_rpm$RPM_fb,miR_merge_rpm$RPM_leaf)$p.value
ks.test(miR_merge_rpm$RPM_fb,miR_merge_rpm$RPM_d234)$p.value

ks.test(tasiR_merge_rpm$RPM_fb,tasiR_merge_rpm$RPM_leaf)$p.value
ks.test(tasiR_merge_rpm$RPM_fb,tasiR_merge_rpm$RPM_d234)$p.value

ks.test(siR_21_merge_rpm$RPM_fb,siR_21_merge_rpm$RPM_leaf)$p.value
ks.test(siR_21_merge_rpm$RPM_fb,siR_21_merge_rpm$RPM_d234)$p.value

ks.test(siR_24_merge_rpm$RPM_fb,siR_24_merge_rpm$RPM_leaf)$p.value
ks.test(siR_24_merge_rpm$RPM_fb,siR_24_merge_rpm$RPM_d234)$p.value

#plot mpu-based violin plots

pdf("figure.2/sRNA_violin_mpu.pdf", family='Arial', useDingbats=F)
par(mfrow=c(2,2), las=1)
vioplot(miR_merge_mpu_log$MPU_leaf,miR_merge_mpu_log$MPU_fb,miR_merge_mpu_log$MPU_d234, ylim=c(4,7.5), col='white', names=name_v)
vioplot(miR_merge_mpu_log$MPU_leaf,col=col_v[1], at=1, add=T)
vioplot(miR_merge_mpu_log$MPU_fb,col=col_v[2], at=2, add=T)
vioplot(miR_merge_mpu_log$MPU_d234,col=col_v[3], at=3, add=T)
text(0.4,7.4,labels=c("miRNAs"), pos=4)
abline(h=median(miR_merge_mpu_log$MPU_fb), lty=2, lwd=1, col="grey")
#tasiRNAs
vioplot(tasiR_merge_mpu_log$MPU_leaf,tasiR_merge_mpu_log$MPU_fb,tasiR_merge_mpu_log$MPU_d234, ylim=c(4,7.5), col='white', names=name_v)
vioplot(tasiR_merge_mpu_log$MPU_leaf,col=col_v[1], at=1, add=T)
vioplot(tasiR_merge_mpu_log$MPU_fb,col=col_v[2], at=2, add=T)
vioplot(tasiR_merge_mpu_log$MPU_d234,col=col_v[3], at=3, add=T)
text(0.4,7.4,labels=c("tasiRNAs"), pos=4)
abline(h=median(tasiR_merge_mpu_log$MPU_fb), lty=2, lwd=1, col="grey")
#20-22 nt siRNAs
vioplot(siR_21_merge_mpu_log$MPU_leaf,siR_21_merge_mpu_log$MPU_fb,siR_21_merge_mpu_log$MPU_d234, ylim=c(4,7.5), col='white', names=name_v)
vioplot(siR_21_merge_mpu_log$MPU_leaf,col=col_v[1], at=1, add=T)
vioplot(siR_21_merge_mpu_log$MPU_fb,col=col_v[2], at=2, add=T)
vioplot(siR_21_merge_mpu_log$MPU_d234,col=col_v[3], at=3, add=T)
text(0.4,7.4,labels=c("20-22 nt siRNAs"), pos=4)
abline(h=median(siR_21_merge_mpu_log$MPU_fb), lty=2, lwd=1, col="grey")
#23-24 nt siRNAs
vioplot(siR_24_merge_mpu_log$MPU_leaf,siR_24_merge_mpu_log$MPU_fb,siR_24_merge_mpu_log$MPU_d234, ylim=c(4,7.5), col='white', names=name_v)
vioplot(siR_24_merge_mpu_log$MPU_leaf,col=col_v[1], at=1, add=T)
vioplot(siR_24_merge_mpu_log$MPU_fb,col=col_v[2], at=2, add=T)
vioplot(siR_24_merge_mpu_log$MPU_d234,col=col_v[3], at=3, add=T)
text(0.4,7.4,labels=c("23-24 nt siRNAs"), pos=4)
abline(h=median(siR_24_merge_mpu_log$MPU_fb), lty=2, lwd=1, col="grey")
dev.off()

ks.test(miR_merge_mpu$MPU_fb,miR_merge_mpu$MPU_leaf)$p.value
ks.test(miR_merge_mpu$MPU_fb,miR_merge_mpu$MPU_d234)$p.value

ks.test(tasiR_merge_mpu$MPU_fb,tasiR_merge_mpu$MPU_leaf)$p.value
ks.test(tasiR_merge_mpu$MPU_fb,tasiR_merge_mpu$MPU_d234)$p.value

ks.test(siR_21_merge_mpu$MPU_fb,siR_21_merge_mpu$MPU_leaf)$p.value
ks.test(siR_21_merge_mpu$MPU_fb,siR_21_merge_mpu$MPU_d234)$p.value

ks.test(siR_24_merge_mpu$MPU_fb,siR_24_merge_mpu$MPU_leaf)$p.value
ks.test(siR_24_merge_mpu$MPU_fb,siR_24_merge_mpu$MPU_d234)$p.value


#scatterplots
##determine how many up and down regulated miRNAs (FDR < 0.05 and >= 2-fold difference) 
fb_vs_leaf = getVals_sig('col0_fb','col0_leaf')
fb_vs_d234 = getVals_sig('col0_fb','d234_fb')
###extract individual list components
fb_vs_leaf_rpm_up = log(fb_vs_leaf[[1]],10)
fb_vs_leaf_rpm_down = log(fb_vs_leaf[[2]],10)
fb_vs_leaf_mpu_up = fb_vs_leaf[[3]]
fb_vs_leaf_mpu_down = fb_vs_leaf[[4]]

fb_vs_d234_rpm_up = log(fb_vs_d234[[1]],10)
fb_vs_d234_rpm_down = log(fb_vs_d234[[2]],10)
fb_vs_d234_mpu_up = fb_vs_d234[[3]]
fb_vs_d234_mpu_down = fb_vs_d234[[4]]

#plot rpm-based scatterplots and add positively and negatively regulated miRNAs
pdf("figure.2/miRNA_scatter_leaf_rpm.pdf", family='Arial', useDingbats=F)
plot(miR_merge_rpm_log$RPM_fb,miR_merge_rpm_log$RPM_leaf, xlim= c(0,4.7), ylim=c(0,4.7), col='black', xlab="Col-0 flowers RPM (log10)", ylab="Col-0 leaves RPM (log10)", las=1)
points(fb_vs_leaf_rpm_up$col0_fb,fb_vs_leaf_rpm_up$col0_leaf, col="red", pch=19)
points(fb_vs_leaf_rpm_down$col0_fb,fb_vs_leaf_rpm_down$col0_leaf, col="blue", pch=19)
abline(a = 0, b = 1, lty=2, lwd=1, col="grey")
dev.off()

pdf("figure.2/miRNA_scatter_d234_rpm.pdf", family='Arial', useDingbats=F)
plot(miR_merge_rpm_log$RPM_fb,miR_merge_rpm_log$RPM_d234, xlim= c(0,4.7), ylim=c(0,4.7), col='black', xlab="Col-0 flowers RPM (log10)", ylab="dcl234 flowers RPM (log10)", las=1)
points(fb_vs_d234_rpm_up$col0_fb,fb_vs_d234_rpm_up$d234_fb, col="red", pch=19)
points(fb_vs_d234_rpm_down$col0_fb,fb_vs_d234_rpm_down$d234_fb, col="blue", pch=19)
abline(a = 0, b = 1, lty=2, lwd=1, col="grey")
dev.off()

#plot mpu-based scatterplots and add positively and negatively regulated miRNAs
pdf("figure.2/miRNA_scatter_leaf_mpu.pdf", family='Arial', useDingbats=F)
plot(miR_merge_mpu_log$MPU_fb,miR_merge_mpu_log$MPU_leaf, xlim= c(4,7.2), ylim=c(4,7.2), col='black', xlab="Col-0 flowers MPU (log10)", ylab="Col-0 leaves MPU (log10)", las=1)
points(fb_vs_leaf_mpu_up$col0_fb,fb_vs_leaf_mpu_up$col0_leaf, col="red", pch=19)
points(fb_vs_leaf_mpu_down$col0_fb,fb_vs_leaf_mpu_down$col0_leaf, col="blue", pch=19)
abline(a = 0, b = 1, lty=2, lwd=1, col="grey")
dev.off()

pdf("figure.2/miRNA_scatter_d234_mpu.pdf", family='Arial', useDingbats=F)
plot(miR_merge_mpu_log$MPU_fb,miR_merge_mpu_log$MPU_d234, xlim= c(4.5,7.5), ylim=c(4.5,7.5), col='black', xlab="Col-0 flowers MPU (log10)", ylab="dcl234 flowers MPU (log10)", las=1)
points(fb_vs_d234_mpu_up$col0_fb,fb_vs_d234_mpu_up$d234_fb, col="red", pch=19)
points(fb_vs_d234_mpu_down$col0_fb,fb_vs_d234_mpu_down$d234_fb, col="blue", pch=19)
abline(a = 0, b = 1, lty=2, lwd=1, col="grey")
dev.off()

