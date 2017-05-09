#generate.rpm.and.mpu.strips.R
#plots all strip charts of miRNAs, tasiRNAs, 20-22 nt siRNAs and 23-24 nt siRNAs as shown in Supplementary Figure 2D-2E

#install.packages("extrafont")
library("extrafont")
#font_import()
loadfonts()

#install.packages("RColorBrewer")
library("RColorBrewer")

data = read.table("methods.data.analysis/small.RNA.Seq/data.for.graphs/total.rpm.mlc.table", header=TRUE, sep="\t", row.names=1, strip.white=T)

#extract subsets for individual plotting
##RPM subsets; note: divide by 1000 for more meaningful values
fb_rpm_mir = c(data$miR_rpm[1],data$miR_rpm[2])/1000
leaf_rpm_mir = c(data$miR_rpm[3],data$miR_rpm[4])/1000
d234_rpm_mir = c(data$miR_rpm[5],data$miR_rpm[6])/1000

fb_rpm_sir24 = c(data$siR_24_rpm[1],data$siR_24_rpm[2])/1000
leaf_rpm_sir24 = c(data$siR_24_rpm[3],data$siR_24_rpm[4])/1000
d234_rpm_sir24 = c(data$siR_24_rpm[5],data$siR_24_rpm[6])/1000

fb_rpm_sir21 = c(data$siR_21_rpm[1],data$siR_21_rpm[2])/1000
leaf_rpm_sir21 = c(data$siR_21_rpm[3],data$siR_21_rpm[4])/1000
d234_rpm_sir21 = c(data$siR_21_rpm[5],data$siR_21_rpm[6])/1000

fb_rpm_tas = c(data$tasiR_rpm[1],data$tasiR_rpm[2])/1000
leaf_rpm_tas = c(data$tasiR_rpm[3],data$tasiR_rpm[4])/1000
d234_rpm_tas = c(data$tasiR_rpm[5],data$tasiR_rpm[6])/1000

##MPU subsets; note: divide by one million for more meaningful values
fb_mlc_mir = c(data$miR_mlc[1],data$miR_mlc[2])/1000000
leaf_mlc_mir = c(data$miR_mlc[3],data$miR_mlc[4])/1000000
d234_mlc_mir = c(data$miR_mlc[5],data$miR_mlc[6])/1000000

fb_mlc_sir24 = c(data$siR_24_mlc[1],data$siR_24_mlc[2])/1000000
leaf_mlc_sir24 = c(data$siR_24_mlc[3],data$siR_24_mlc[4])/1000000
d234_mlc_sir24 = c(data$siR_24_mlc[5],data$siR_24_mlc[6])/1000000

fb_mlc_sir21 = c(data$siR_21_mlc[1],data$siR_21_mlc[2])/1000000
leaf_mlc_sir21 = c(data$siR_21_mlc[3],data$siR_21_mlc[4])/1000000
d234_mlc_sir21 = c(data$siR_21_mlc[5],data$siR_21_mlc[6])/1000000

fb_mlc_tas = c(data$tasiR_mlc[1],data$tasiR_mlc[2])/1000000
leaf_mlc_tas = c(data$tasiR_mlc[3],data$tasiR_mlc[4])/1000000
d234_mlc_tas = c(data$tasiR_mlc[5],data$tasiR_mlc[6])/1000000

col_v <- brewer.pal(4,"Paired") 
col_v = c(col_v[4],col_v[3],col_v[1])

pdf("supplemental.figure.2/sRNA.total.rpm.strip.pdf", family='Arial', useDingbats=F)
stripchart(list(leaf_rpm_mir,fb_rpm_mir,d234_rpm_mir,leaf_rpm_tas,fb_rpm_tas,d234_rpm_tas,leaf_rpm_sir21,fb_rpm_sir21,d234_rpm_sir21,leaf_rpm_sir24,fb_rpm_sir24,d234_rpm_sir24), pch=20, las=1, ylab="Thousands of reads per million genome-matching reads", vertical=T, method='jitter', col='white')
stripchart(leaf_rpm_mir,col=col_v[1], at=1, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_rpm_mir,col=col_v[2], at=2, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_rpm_mir,col=col_v[3], at=3, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_rpm_tas,col=col_v[1], at=4, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_rpm_tas,col=col_v[2], at=5, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_rpm_tas,col=col_v[3], at=6, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_rpm_sir21,col=col_v[1], at=7, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_rpm_sir21,col=col_v[2], at=8, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_rpm_sir21,col=col_v[3], at=9, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_rpm_sir24,col=col_v[1], at=10, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_rpm_sir24,col=col_v[2], at=11, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_rpm_sir24,col=col_v[3], at=12, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
legend_coords = xy.coords(1,200)
legend(legend_coords,c('Col-0 leaves','Col-0 flowers','dcl234 flowers'), pch=19, col=col_v, pt.cex=2, bty="n")
dev.off()

t.test(leaf_rpm_mir,fb_rpm_mir)$p.value
t.test(d234_rpm_mir,fb_rpm_mir)$p.value

t.test(leaf_rpm_tas,fb_rpm_tas)$p.value
t.test(d234_rpm_tas,fb_rpm_tas)$p.value

t.test(leaf_rpm_sir21,fb_rpm_sir21)$p.value
t.test(d234_rpm_sir21,fb_rpm_sir21)$p.value

t.test(leaf_rpm_sir24,fb_rpm_sir24)$p.value
t.test(d234_rpm_sir24,fb_rpm_sir24)$p.value


pdf("supplemental.figure.2/sRNA.total.mlc.strip.pdf", family='Arial', useDingbats=F)
stripchart(list(leaf_mlc_mir,fb_rpm_mir,d234_mlc_mir,leaf_mlc_tas,fb_mlc_tas,d234_mlc_tas,leaf_mlc_sir21,fb_mlc_sir21,d234_mlc_sir21,leaf_mlc_sir24,fb_mlc_sir24,d234_mlc_sir24), pch=20, las=1, ylab="Millions of molecules per ug total RNA", vertical=T, method='jitter', col='white')
stripchart(leaf_mlc_mir,col=col_v[1], at=1, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_mlc_mir,col=col_v[2], at=2, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_mlc_mir,col=col_v[3], at=3, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_mlc_tas,col=col_v[1], at=4, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_mlc_tas,col=col_v[2], at=5, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_mlc_tas,col=col_v[3], at=6, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_mlc_sir21,col=col_v[1], at=7, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_mlc_sir21,col=col_v[2], at=8, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_mlc_sir21,col=col_v[3], at=9, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(leaf_mlc_sir24,col=col_v[1], at=10, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(fb_mlc_sir24,col=col_v[2], at=11, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
stripchart(d234_mlc_sir24,col=col_v[3], at=12, add=T, pch=20, cex=2, vertical=T, method='jitter', jitter=0.3)
legend_coords = xy.coords(1,80)
legend(legend_coords,c('Col-0 leaves','Col-0 flowers','dcl234 flowers'), pch=19, col=col_v, pt.cex=2, bty="n")
dev.off()

t.test(leaf_mlc_mir,fb_mlc_mir)$p.value
t.test(d234_mlc_mir,fb_mlc_mir)$p.value

t.test(leaf_mlc_tas,fb_mlc_tas)$p.value
t.test(d234_mlc_tas,fb_mlc_tas)$p.value

t.test(leaf_mlc_sir21,fb_mlc_sir21)$p.value
t.test(d234_mlc_sir21,fb_mlc_sir21)$p.value

t.test(leaf_mlc_sir24,fb_mlc_sir24)$p.value
t.test(d234_mlc_sir24,fb_mlc_sir24)$p.value
