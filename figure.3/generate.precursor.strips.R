#generate.precursor.strips.R
#creates strip charts of mature/precursor ratios for miRNAs and tasiRNAs as shown in Fig. 3B

#install.packages("extrafont")
library("extrafont")
#font_import()
loadfonts()

#install.packages("RColorBrewer")
library("RColorBrewer")

setwd('/Volumes/nodine/lab/members/michael/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/')  #user will have to change this to their working directory

master = read.table("methods.data.analysis/master_genes_Araport.txt", header=TRUE, sep='\t', row.names=1, strip.white=T)

Na = 6.02214179 * 10e23 * 1e-18

getVals_prec <- function(sample,type,mix,scaling_factor){
  data = read.table(paste("methods.data.analysis/small.RNA.Seq/data.for.graphs/",sample,'_',type,'_prec',sep=''), header=T, sep='\t', row.names=1, strip.white=T)
  
  #select mature miRNAs that are expressed; i.e. that have an RPM >= 1
  data_sub = subset(data, sRNA_rpm >= 1)
  
  #now extract TPM values for AGIs
  ##take subset of master corresponding to sample
  master_sub = subset(master, select = sample)  
  prec = merge(master_sub,data_sub, by='row.names')
  rownames(prec) = prec$Row.names
  prec_exp = subset(prec, prec[,sample] >= 1, select = 2:3)
  colnames(prec_exp) = c("TPM","sRNA_rpm")

  spike = read.table("methods.data.analysis/Analysis.txt", header=TRUE, sep="\t", row.names=2, strip.white=T)
  spike = subset(spike, select=paste("concentration.in.Mix.",mix,"..attomoles.ul.",sep=""))
  colnames(spike) = "molecules"
  
  #multiply by Avagadro's number, scaling factor and 1/1000 to determine number of molecules for each spike-in added
  #scaling_factor was the same (0.0002) all samples used in this study; scaling_factor = (ul of ERCC mix added) * (1/dilution of total ug that spike-in was added to) * (1/ERCC dilution) * (1/total ug)
  # = (1) * (1/50) * (1/200) * (1/0.5)
  
  spike = spike * Na * scaling_factor
  
  #select ERCCs with values >= 1 TPM
  ERCC_all = merge(master_sub,spike, by = "row.names")
  row.names(ERCC_all) = ERCC_all$Row.names
  ERCC_exp = subset(ERCC_all, ERCC_all[,sample] >= 1, select = 2:3)
  colnames(ERCC_exp) = c("TPM","molecules")

  fit_log = lm(log(ERCC_exp$molecules,10) ~ log(ERCC_exp$TPM,10))
  
  #now get sRNA molecules  
  data_sRNA = read.table(paste("methods.data.analysis/small.RNA.Seq/",sample,"/smRNA_spikeIns/doseResponseTable_noTransform",sep=""), header=TRUE, sep="\t", row.names=1, strip.white=T)

  #now fit linear model and use to predict mpe
  fit_log_sRNA = lm(log(data_sRNA$mlcs,10) ~ log(data_sRNA$rpm..mean.,10))
    
  #now apply model to precursor TPM and sRNA rpm values, and estimate individual molecules  
  coeffs_mRNA = coefficients(fit_log)
  coeffs_sRNA = coefficients(fit_log_sRNA)
  
  absolute = cbind(prec_exp$TPM, 10^(coeffs_mRNA[1] + coeffs_mRNA[2]*log(prec_exp$TPM,10)), prec_exp$sRNA_rpm, 10^(coeffs_sRNA[1] + coeffs_sRNA[2]*log(prec_exp$sRNA_rpm,10)))   #, mean(c(coeffs_1[1] + coeffs_1[2]*sub$biorep1,coeffs_2[1] + coeffs_2[2]*sub$biorep2)))
  dimnames(absolute)=list(c(row.names(prec_exp)),c("TPM","mRNA_mlcs","RPM","sRNA_mlcs"))
    
  return(absolute)
}

col0_fb1_tas = getVals_prec('col0_fb1','tas', 1, 0.0002)
col0_fb2_tas = getVals_prec('col0_fb2','tas', 1, 0.0002)

d234_fb1_tas = getVals_prec('d234_fb1','tas', 2, 0.0002)
d234_fb2_tas = getVals_prec('d234_fb2','tas', 2, 0.0002)

col0_leaf1_tas = getVals_prec('col0_leaf1','tas', 1, 0.0002)
col0_leaf2_tas = getVals_prec('col0_leaf2','tas', 1, 0.0002)

col0_fb1_mir = getVals_prec('col0_fb1','miR', 1, 0.0002)
col0_fb2_mir = getVals_prec('col0_fb2','miR', 1, 0.0002)

d234_fb1_mir = getVals_prec('d234_fb1','miR', 2, 0.0002)
d234_fb2_mir = getVals_prec('d234_fb2','miR', 2, 0.0002)

col0_leaf1_mir = getVals_prec('col0_leaf1','miR', 1, 0.0002)
col0_leaf2_mir = getVals_prec('col0_leaf2','miR', 1, 0.0002)

##col0_leaf_mir
col0_leaf_mir = merge(col0_leaf1_mir,col0_leaf2_mir, by = "row.names")
row.names(col0_leaf_mir) = col0_leaf_mir$Row.names
col0_leaf_mir = subset(col0_leaf_mir, select = c(3,5,7,9))
colnames(col0_leaf_mir) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")
##col0_fb_mir
col0_fb_mir = merge(col0_fb1_mir,col0_fb2_mir, by = "row.names")
row.names(col0_fb_mir) = col0_fb_mir$Row.names
col0_fb_mir = subset(col0_fb_mir, select = c(3,5,7,9))
colnames(col0_fb_mir) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")
##d234_fb_mir
d234_fb_mir = merge(d234_fb1_mir,d234_fb2_mir, by = "row.names")
row.names(d234_fb_mir) = d234_fb_mir$Row.names
d234_fb_mir = subset(d234_fb_mir, select = c(3,5,7,9))
colnames(d234_fb_mir) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")

##col0_leaf_tas
col0_leaf_tas = merge(col0_leaf1_tas,col0_leaf2_tas, by = "row.names")
row.names(col0_leaf_tas) = col0_leaf_tas$Row.names
col0_leaf_tas = subset(col0_leaf_tas, select = c(3,5,7,9))
colnames(col0_leaf_tas) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")
##col0_fb_tas
col0_fb_tas = merge(col0_fb1_tas,col0_fb2_tas, by = "row.names")
row.names(col0_fb_tas) = col0_fb_tas$Row.names
col0_fb_tas = subset(col0_fb_tas, select = c(3,5,7,9))
colnames(col0_fb_tas) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")
##d234_fb_tas
d234_fb_tas = merge(d234_fb1_tas,d234_fb2_tas, by = "row.names")
row.names(d234_fb_tas) = d234_fb_tas$Row.names
d234_fb_tas = subset(d234_fb_tas, select = c(3,5,7,9))
colnames(d234_fb_tas) = c("mRNA_mlcs_rep1","sRNA_mlcs_rep1","mRNA_mlcs_rep2","sRNA_mlcs_rep2")

col_v <- brewer.pal(4,"Paired") 
col_v = c(col_v[4],col_v[3],col_v[1])
name_v = c("leaves","flowers","d234","","leaves","flowers","d234")

col0_leaf_mir_vals = (log(col0_leaf_mir[,'sRNA_mlcs_rep1']/col0_leaf_mir[,'mRNA_mlcs_rep1'],2) + log(col0_leaf_mir[,'sRNA_mlcs_rep2']/col0_leaf_mir[,'mRNA_mlcs_rep2'],2))/2
col0_fb_mir_vals = (log(col0_fb_mir[,'sRNA_mlcs_rep1']/col0_fb_mir[,'mRNA_mlcs_rep1'],2) + log(col0_fb_mir[,'sRNA_mlcs_rep2']/col0_fb_mir[,'mRNA_mlcs_rep2'],2))/2
d234_fb_mir_vals = (log(d234_fb_mir[,'sRNA_mlcs_rep1']/d234_fb_mir[,'mRNA_mlcs_rep1'],2) + log(d234_fb_mir[,'sRNA_mlcs_rep2']/d234_fb_mir[,'mRNA_mlcs_rep2'],2))/2

col0_leaf_tas_vals = (log(col0_leaf_tas[,'sRNA_mlcs_rep1']/col0_leaf_tas[,'mRNA_mlcs_rep1'],2) + log(col0_leaf_tas[,'sRNA_mlcs_rep2']/col0_leaf_tas[,'mRNA_mlcs_rep2'],2))/2
col0_fb_tas_vals = (log(col0_fb_tas[,'sRNA_mlcs_rep1']/col0_fb_tas[,'mRNA_mlcs_rep1'],2) + log(col0_fb_tas[,'sRNA_mlcs_rep2']/col0_fb_tas[,'mRNA_mlcs_rep2'],2))/2
d234_fb_tas_vals = (log(d234_fb_tas[,'sRNA_mlcs_rep1']/d234_fb_tas[,'mRNA_mlcs_rep1'],2) + log(d234_fb_tas[,'sRNA_mlcs_rep2']/d234_fb_tas[,'mRNA_mlcs_rep2'],2))/2

pdf("figure.3/sRNA.precursor.strip.pdf", family='Arial', useDingbats=F)
stripchart(list(col0_leaf_mir_vals,col0_fb_mir_vals,d234_fb_mir_vals,c(0),col0_leaf_tas_vals,col0_fb_tas_vals,d234_fb_tas_vals), pch=20, las=1, ylab="Mature/precursor (log2)", vertical=T, method='jitter', col='white', group.names=name_v)
stripchart(col0_leaf_mir_vals,col=col_v[1], at=1, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(col0_fb_mir_vals,col=col_v[2], at=2, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_fb_mir_vals,col=col_v[3], at=3, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(col0_leaf_tas_vals,col=col_v[1], at=5, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(col0_fb_tas_vals,col=col_v[2], at=6, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_fb_tas_vals,col=col_v[3], at=7, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
arrows(1 - 0.25,median(col0_leaf_mir_vals),1 + 0.25,median(col0_leaf_mir_vals),length=0,lwd=2,col="black")
arrows(2 - 0.25,median(col0_fb_mir_vals),2 + 0.25,median(col0_fb_mir_vals),length=0,lwd=2,col="black")
arrows(3 - 0.25,median(d234_fb_mir_vals),3 + 0.25,median(d234_fb_mir_vals),length=0,lwd=2,col="black")
arrows(5 - 0.25,median(col0_leaf_tas_vals),5 + 0.25,median(col0_leaf_tas_vals),length=0,lwd=2,col="black")
arrows(6 - 0.25,median(col0_fb_tas_vals),6 + 0.25,median(col0_fb_tas_vals),length=0,lwd=2,col="black")
arrows(7 - 0.25,median(d234_fb_tas_vals),7 + 0.25,median(d234_fb_tas_vals),length=0,lwd=2,col="black")
legend_coords = xy.coords(5.5,6.8)
legend(legend_coords,c('Col-0 leaves','Col-0 flowers','dcl234 flowers'), pch=15, col=col_v, pt.cex=2.5, bty="n")
dev.off()


ks.test(col0_fb_mir_vals,col0_leaf_mir_vals)$p.value
ks.test(col0_fb_mir_vals,d234_fb_mir_vals)$p.value

ks.test(col0_fb_tas_vals,col0_leaf_tas_vals)$p.value
ks.test(col0_fb_tas_vals,d234_fb_tas_vals)$p.value
