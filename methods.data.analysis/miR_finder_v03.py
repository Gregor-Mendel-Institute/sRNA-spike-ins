#miR_finder_v03.py
#organizes tab-delimited file of sRNA-Seq alignments
#input:
#   1. sample
#   2. total_reads = number of genome-matching reads to use for normalization
#   3. filePath = filepath where files associated with given sample are located
#output:
#   1. "summary" file containing information about output contents and general statistics
#   2. "abundantFound" tab-delimited file of most abundant tags for each MIRNA gene
#   3. "abundantNotMatFound" tab-delimited file of most abundant tag which is not annotated
#   4. "families_byRpm" ordered table of miRNA families found and their RPM levels
#   5. "mature_byRpm" ordered table of mature miRNAs found and their RPM levels
#   6. "matureFound" tab-delimited file of tags contained within each annotated mature miRNA

import sys, os, basicProcessModule_slim, miRNA_module_slim, operator

sample = sys.argv[1].strip()
total_reads = float(sys.argv[2].strip())
filePath = sys.argv[3].strip()

dataPath = filePath + '/' + sample + '/'
gmr = total_reads

outPath = dataPath + '/miRNA/'

cmd = 'mkdir -p ' + outPath
os.system(cmd)


#make ghc of reads mapping to sense strand of MIRNA genes

print 'Making ghc...'

ghc = basicProcessModule_slim.makeGHC_fromHitsFile_v05(dataPath,'mi','all','sense')

#make MIRNA locus collection to search for overlaps

lc_mi = basicProcessModule_slim.makeLocusCollection('mi')

sumFile = outPath + 'summary'
sfh = open(sumFile,'w')

print >>sfh,'The following files are contained within this directory:'
print >>sfh,'...matureFound: contains values for reads contained within mature miRNAs (+/- 2 bp) annotated in TAIR10/miRBase21'
print >>sfh,'...abundantFound: contains values for most abundant tag overlapping each annotated MIRNA locus'
print >>sfh,'...abundantNotMatureFound: contains values for most abundant tag overlapping each annotated MIRNA locus if not contained in annotated mature'
print >>sfh,'...mature_byRpm: contains ordered list of mature miRNAs'
print >>sfh,'...families_byRpm: contains ordered list of families'

print 'Number of annotated MIRNA loci:',len(lc_mi.getLoci())
print >>sfh,'Number of annotated MIRNA loci:',len(lc_mi.getLoci())

print 'Number of alignments to sense strand of MIRNA genes:',len(ghc.getHits())
print >>sfh,'Number of alignments to sense strand of MIRNA genes:',len(ghc.getHits())

print 'Number of (hit-normalized) reads (rpm) that map to sense strand of MIRNA genes:',round(ghc.get_readNum_normHit(),2),'(' + str(round(ghc.get_readNum_normHit()/gmr,2)) + ')'
print >>sfh,'Number of (hit-normalized) reads (rpm) that map to sense strand of MIRNA genes:',round(ghc.get_readNum_normHit(),2),'(' + str(round(ghc.get_readNum_normHit()/gmr,2)) + ')'


#search for overlap

hitDict = miRNA_module_slim.getOverlappingLoci_miRNA(ghc,lc_mi,'sense')

print 'Number of MIRNA genes with sRNA reads:',len(hitDict.keys())
print >>sfh,'Number of MIRNA genes with sRNA reads:',len(hitDict.keys())

matDict = miRNA_module_slim.getContained_mature(hitDict)

abundDict = miRNA_module_slim.getAbund(hitDict)

print 'Number of mature miRNAs with corresponding reads:',len(matDict.keys())
print 'Number of reads determined to be most abundant for given MIRNA:',len(abundDict.keys())
print >>sfh,'Number of mature miRNAs with corresponding reads:',len(matDict.keys())
print >>sfh,'Number of reads determined to be most abundant for given MIRNA:',len(abundDict.keys())


matFile = outPath + 'matureFound'
abundFile = outPath + 'abundantFound'
abundNotMatFile = outPath + 'abundantNotMatFound'
matFile_ordered = outPath + 'mature_byRpm'
famFile = outPath + 'families_byRpm'

matDict2 = {}
famDict = {}

annoFamDict = miRNA_module_slim.getAnnoFamDict()

mfh = open(matFile,'w')
afh = open(abundFile,'w')
nfh = open(abundNotMatFile,'w')
mfh2 = open(matFile_ordered,'w')
ffh = open(famFile,'w')
	
print >>mfh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'
print >>afh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'
print >>nfh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'

for locus in matDict.keys():
    mat = locus.geneID()      #locus.short_desc()
    fam = annoFamDict[mat]  #locus.long_desc()
    count = 0
    rawCount = 0
    ghc_mat = matDict[locus]
    for hit in ghc_mat.getHits():
        count = count + hit.readNum_normRpm()
        rawCount = rawCount  + hit.readNum_normHit()
	print >>mfh,mat,'\t',hit.chromo(),'\t',hit.start(),'\t',hit.end(),'\t',hit.strand(),'\t',hit.seq(),'\t',hit.readNum_raw(),'\t',hit.readNum_normHit(),'\t',hit.readNum_normRpm()

    if len(ghc_mat.getHits()) > 0:
        if mat in matDict2.keys():
            matCount = matDict2[mat] + count
            matDict2[mat] = (matCount)
        else:
            matDict2[mat] = (count)

    if len(ghc_mat.getHits()) > 0:
        if fam in famDict.keys():
            famCount = famDict[fam][0] + count
	    rawCount2 = famDict[fam][1] + rawCount
            famDict[fam] = (famCount,rawCount2)
        else:
            famDict[fam] = (count,rawCount)

matDict_sorted = sorted(matDict2.items(),key=operator.itemgetter(1), reverse=True)
famDict_sorted = sorted(famDict.items(),key=operator.itemgetter(1), reverse=True)

miRNA_total = 0

for fam in famDict.keys():
    miRNA_total = miRNA_total + famDict[fam][1]  

print "Number of reads contained within miRBase21/TAIR10 annotated miRNAs:",miRNA_total

print >>mfh2,'mature','\t','rpm'

#print >>ffh,"####Number of reads conained within miRBase21/TAIR10 annotated miRNAs:",miRNA_total
print >>ffh,'family','\t','rpm','\t','rpmir'

for mat in matDict_sorted:
    print >>mfh2,mat[0],'\t',mat[1]

for fam in famDict_sorted:
    print >>ffh,fam[0],'\t',fam[1][0],'\t',float(fam[1][1])/miRNA_total
        
for locus in abundDict.keys():
    maxHit = abundDict[locus]
    print >>afh,locus.common(),'\t',maxHit.chromo(),'\t',maxHit.start(),'\t',maxHit.end(),'\t',maxHit.strand(),'\t',maxHit.seq(),'\t',maxHit.readNum_raw(),'\t',maxHit.readNum_normHit(),'\t',maxHit.readNum_normRpm()


hitList = []

for locus in matDict.keys():
    ghc_mat = matDict[locus]
    for hit in ghc_mat.getHits():
        hit = hit.chromo() + '_' + str(hit.start()) + '_' + str(hit.end()) + '_' + hit.strand()
        hitList.append(hit)

notInCount = 0

for locus in abundDict.keys():
    maxHit = abundDict[locus]
    hit = maxHit.chromo() + '_' + str(maxHit.start()) + '_' + str(maxHit.end()) + '_' + maxHit.strand()
    if hit not in hitList:
        notInCount += 1
        print >>nfh,locus.common(),'\t',maxHit.chromo(),'\t',maxHit.start(),'\t',maxHit.end(),'\t',maxHit.strand(),'\t',maxHit.seq(),'\t',maxHit.readNum_raw(),'\t',maxHit.readNum_normHit(),'\t',maxHit.readNum_normRpm()

print 'Number of MIRNA loci with an unannotated tag being the most abundant:',notInCount,round(notInCount/float(len(hitDict.keys())),2)
print >>sfh,'Number of MIRNA loci with an unannotated tag being the most abundant:',notInCount,round(notInCount/float(len(hitDict.keys())),2)

mfh.close()
mfh2.close()
afh.close()
nfh.close()
ffh.close()
sfh.close()
