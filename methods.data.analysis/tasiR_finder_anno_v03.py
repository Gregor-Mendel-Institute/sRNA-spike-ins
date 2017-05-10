#tasiR_finder_anno_v03.py
#organizes tab-delimited file of sRNA-Seq alignments
#input:
#   1. sample
#   2. total_reads = number of genome-matching reads to use for normalization
#   3. filePath = filepath where files associated with given sample are located
#output:
#   1. "summary" file containing information about output contents and general statistics
#   2. "abundantFound" tab-delimited file of most abundant tags for each tasiRNA gene
#   3. "abundantNotMatFound" tab-delimited file of most abundant tag which is not annotated
#   4. "families_byRpm" ordered table of tasiRNA families found and their RPM levels
#   5. "mature_byRpm" ordered table of mature tasiRNAs found and their RPM levels
#   6. "matureFound" tab-delimited file of tags contained within each annotated mature tasiRNA



import sys, os, basicProcessModule_slim, miRNA_module_slim, operator

sample = sys.argv[1].strip()
total_reads = float(sys.argv[2].strip())
filePath = sys.argv[3].strip()

dataPath = filePath + '/' + sample + '/'
gmr = total_reads

outPath = dataPath + '/tasiRNA/'

cmd = 'mkdir -p ' + outPath
os.system(cmd)

#make ghc of reads mapping to either strand of TASIRNA genes

print 'Making ghc...'

ghc = basicProcessModule_slim.makeGHC_fromHitsFile_v05(dataPath,'tas','all','both')

#make TASIRNA locus collection to search for overlaps

lc_tas = basicProcessModule_slim.makeLocusCollection('tas')

sumFile = outPath + 'summary'
sfh = open(sumFile,'w')

print >>sfh,'The following files are contained within this directory:'
print >>sfh,'...matureFound: contains values for reads contained within tasiRNAs (+/- 2 bp) annotated by Allen et al. (2005) Cell and tasiRNAdb'
print >>sfh,'...abundantFound: contains values for most abundant tag overlapping each annotated TASIRNA locus'
print >>sfh,'...abundantNotMatureFound: contains values for most abundant tag overlapping each annotated TASIRNA locus if not contained in annotated tasiRNA'
print >>sfh,'...mature_byRpm: contains ordered list of mature tasiRNAs'
print >>sfh,'...families_byRpm: contains ordered list of families'

print 'Number of annotated TASIRNA loci:',len(lc_tas.getLoci())
print >>sfh,'Number of annotated TASIRNA loci:',len(lc_tas.getLoci())

print 'Number of alignments to either strand of TASIRNA genes:',len(ghc.getHits())
print >>sfh,'Number of alignments to either strand of TASIRNA genes:',len(ghc.getHits())

print 'Number of (hit-normalized) reads (rpm) that map to either strand of TASIRNA genes:',round(ghc.get_readNum_normHit(),2),'(' + str(round(ghc.get_readNum_normHit()/gmr,2)) + ')'
print >>sfh,'Number of (hit-normalized) reads (rpm) that map to either strand of TASIRNA genes:',round(ghc.get_readNum_normHit(),2),'(' + str(round(ghc.get_readNum_normHit()/gmr,2)) + ')'

readNumTotal = ghc.get_readNum_normRpm()

hitDict = miRNA_module_slim.getOverlappingLoci_miRNA(ghc,lc_tas,'both')

print 'Number of TASIRNA genes with sRNA reads:',len(hitDict.keys())
print >>sfh,'Number of TASIRNA genes with sRNA reads:',len(hitDict.keys())

matDict = miRNA_module_slim.getContained_tasi(hitDict)

abundDict = miRNA_module_slim.getAbund(hitDict)

print 'Number of tasiRNAs with corresponding reads:',len(matDict.keys())
print 'Number of reads determined to be most abundant for given TASIRNA:',len(abundDict.keys())
print >>sfh,'Number of tasiRNAs with corresponding reads:',len(matDict.keys())
print >>sfh,'Number of reads determined to be most abundant for given TASIRNA:',len(abundDict.keys())

matFile = outPath + 'matureFound'
abundFile = outPath + 'abundantFound'
abundNotMatFile = outPath + 'abundantNotMatFound'
matFile_ordered = outPath + 'mature_byRpm'
famFile = outPath + 'families_byRpm'

matDict2 = {}
famDict = {}

annoFamDict = miRNA_module_slim.getAnnoFamDict_tasi()

#print annoFamDict

mfh = open(matFile,'w')
afh = open(abundFile,'w')
nfh = open(abundNotMatFile,'w')
mfh2 = open(matFile_ordered,'w')
ffh = open(famFile,'w')

print >>mfh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'
print >>afh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'
print >>nfh,'common (miRNA)','\t','chromo','\t','start','\t','end','\t','strand','\t','sequence','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'

annoCount = 0

for locus in matDict.keys():
    count = 0
    mat = locus.geneID()
    fam = annoFamDict[mat.upper()]
    common = locus.short_desc()
    ghc_mat = matDict[locus]
    rawCount = 0
    for hit in ghc_mat.getHits():
        count = count + hit.readNum_normRpm()
        annoCount = annoCount + hit.readNum_normRpm()
        rawCount = rawCount + hit.readNum_normHit()
        print >>mfh,hit.name(),'\t',hit.chromo(),'\t',hit.start(),'\t',hit.end(),'\t',hit.strand(),'\t',hit.seq(),'\t',hit.readNum_raw(),'\t',hit.readNum_normHit(),'\t',hit.readNum_normRpm()

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

tasiRNA_total = 0

for fam in famDict.keys():
    tasiRNA_total = tasiRNA_total + famDict[fam][1]  

print "Number of reads contained within annotated tasiRNAs:",tasiRNA_total

print >>mfh2,'mature','\t','rpm'
print >>ffh,'family','\t','rpm','\t','rptas'

for mat in matDict_sorted:
    print >>mfh2,mat[0],'\t',mat[1]

for fam in famDict_sorted:
    print >>ffh,fam[0],'\t',fam[1][0],'\t',float(fam[1][1])/tasiRNA_total

print 'anno/total tasiRNA rpm:',str(round(annoCount,2)) + '/' + str(round(ghc.get_readNum_normRpm(),2)),round(annoCount/ghc.get_readNum_normRpm(),2)
print >>sfh,'anno/total tasiRNA rpm:',str(round(annoCount,2)) + '/' + str(round(ghc.get_readNum_normRpm(),2)),round(annoCount/ghc.get_readNum_normRpm(),2)

for locus in abundDict.keys():
    maxHit = abundDict[locus]
    print >>afh,locus.short_desc(),'\t',maxHit.chromo(),'\t',maxHit.start(),'\t',maxHit.end(),'\t',maxHit.strand(),'\t',maxHit.seq(),'\t',maxHit.readNum_raw(),'\t',maxHit.readNum_normHit(),'\t',maxHit.readNum_normRpm()

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
    #print hit
    if hit not in hitList:
        notInCount += 1
        print >>nfh,locus.short_desc(),'\t',maxHit.chromo(),'\t',maxHit.start(),'\t',maxHit.end(),'\t',maxHit.strand(),'\t',maxHit.seq(),'\t',maxHit.readNum_raw(),'\t',maxHit.readNum_normHit(),'\t',maxHit.readNum_normRpm()

print 'Number of TASIRNA loci with an unannotated tag being the most abundant:',notInCount,round(notInCount/float(len(hitDict.keys())),2)
print >>sfh,'Number of TASIRNA loci with an unannotated tag being the most abundant:',notInCount,round(notInCount/float(len(hitDict.keys())),2)

mfh.close()
afh.close()
nfh.close()
mfh2.close()
ffh.close()
sfh.close()

