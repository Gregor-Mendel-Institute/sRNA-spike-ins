#make_ghcsAndlDicts.py
#converts allgHits file to individual files containing genomic hit collections (ghc)
#for various gene-types (e.g. miRNAs, tasiRNAs, transposons, etc.)
#input:
#   1. sample
#   2. total_reads = number of genome-matching reads to use for normalization
#   3. mRNA_spikes = whether mRNA spike-ins should be analyzed (yes or no)
#   4. sRNA_spikes = whether sRNA spike-ins should be analyzed (yes or no)
#   5. filePath = filepath where files associated with given sample are located
#output:
#   1. tab-delimited file of ghc for various gene-types (e.g. miRNAs, tasiRNAs, transposons, etc.)
#   2. appends summary information regarding alignment number to various gene-types
#   to existing alignSummary file (output during bwtToAllHits_v03.py)

import sys, basicProcessModule_slim, os

sample = sys.argv[1].strip()
total_reads = float(sys.argv[2].strip())
mRNA_spikes = sys.argv[3].strip()   #'yes' was used in this study
sRNA_spikes = sys.argv[4].strip()   #'yes' was used in this study
filePath = sys.argv[5].strip()


directory = filePath + '/' + sample + '/'

print "Total number of (million) genome-matching reads:",total_reads
        
def printLdict(lDict,name):
    outFile = directory + '/lDicts/' + name
    ofh = open(outFile,'w')

    readNum = 0

    print >>ofh,'locus','\t','tagNum','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm','\t','tagList'

    for locus in lDict.keys():
        if name == 'lDict_mi_sense' or name == 'lDict_mi_antisense':
            if locus.common != '.':
                mir = locus.common().upper()
            else:
                mir = locus.AGI_ID()
            print >>ofh,mir,'\t',lDict[locus][0],'\t',lDict[locus][1],'\t',lDict[locus][2],'\t',lDict[locus][3],'\t',lDict[locus][4]
        else:
            print >>ofh,locus.geneID(),'\t',lDict[locus][0],'\t',lDict[locus][1],'\t',lDict[locus][2],'\t',lDict[locus][3],'\t',lDict[locus][4]            
        readNum = readNum + lDict[locus][2]

    print 'Number of reads in',name,'\t',readNum

def printLdict_spikes(lDict,name):
    outFile = directory + '/lDicts/' + name
    ofh = open(outFile,'w')

    readNum = 0

    print >>ofh,'locus','\t','tagNum','\t','readNum_raw','\t','readNum_normHit','\t','readNum_normRpm'

    for locus in lDict.keys():
        print >>ofh,locus,'\t',lDict[locus][0],'\t',lDict[locus][1],'\t',lDict[locus][2],'\t',lDict[locus][3] 
        readNum = readNum + lDict[locus][2]

    print 'Number of reads in',name,'\t',readNum

def print_gHit(tag,name):
    geneList = tagDict[tag.name()][0]
    typeList = tagDict[tag.name()][1]

    geneList_str = ''
    typeList_str = ''

    for gene in geneList:
        geneList_str = geneList_str + '$' + gene
    for geneType in typeList:
        typeList_str = typeList_str + '$' + geneType

    ofh1 = ofh_dict[name][0]
    ofh2 = ofh_dict[name][1]
    
    print >>ofh1,tag.name(),'\t',tag.seq(),'\t',tag.chromo(),'\t',tag.start(),'\t',tag.end(),'\t',tag.strand(),'\t',tag.readNum_raw(),'\t',tag.readNum_normHit(),'\t',tag.readNum_normRpm(),'\t',geneList_str,'\t',typeList_str

    if len(typeList2) == 1:
        print >>ofh2,tag.name(),'\t',tag.seq(),'\t',tag.chromo(),'\t',tag.start(),'\t',tag.end(),'\t',tag.strand(),'\t',tag.readNum_raw(),'\t',tag.readNum_normHit(),'\t',tag.readNum_normRpm(),'\t',geneList_str,'\t',typeList_str

#first make ghc of all possible tags

print 'Generating ghc_all...'

if mRNA_spikes == 'yes' and sRNA_spikes == 'yes':
    ghc_all,ghc_mRNAspike,ghc_sRNAspike = basicProcessModule_slim.makeGHC_fromAllgHits(directory,mRNA_spikes,sRNA_spikes)
    print 'Total number of genome-matching reads in ghc_all:',ghc_all.get_readNum_normHit()
    print 'Total number of mRNA spike-in matching reads:',ghc_mRNAspike.get_readNum_normHit()
    print 'Total number of sRNA spike-in matching reads:',ghc_sRNAspike.get_readNum_normHit()

#search for overlap with each locus type
##first build lc for all annotated loci

print 'Generating locus collections...'

lc_te = basicProcessModule_slim.makeLocusCollection('te')
lc_mi = basicProcessModule_slim.makeLocusCollection('mi')
lc_tas = basicProcessModule_slim.makeLocusCollection('tas')

print 'Number of annotated te:',len(lc_te.getLoci())
print 'Number of annotated mi:',len(lc_mi.getLoci())
print 'Number of annotated tas:',len(lc_tas.getLoci())

##build agi list for each locus class

print 'Generating AGI lists...'

agi_te = basicProcessModule_slim.makeAGIlist('te')
agi_mi = basicProcessModule_slim.makeAGIlist('mi')
agi_tas = basicProcessModule_slim.makeAGIlist('tas')

##retrieve gHits that overlap each locus class (both, sense, antisense)

###generate tagDict[tag] = ([geneList],[typeList])

print 'Populating tagDict...'

tagDict = {}

for tag in ghc_all.getHits():
    tagDict[tag.name()] = ([],[],[])

print 'Retrieving overlapping gHits...'

print 'te...'
lDict_te_sense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_te,agi_te,'sense',tagDict)
lDict_te_antisense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_te,agi_te,'antisense',tagDict)

print 'mi...'
lDict_mi_sense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_mi,agi_mi,'sense',tagDict)
lDict_mi_antisense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_mi,agi_mi,'antisense',tagDict)

print 'tas...'
lDict_tas_sense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_tas,agi_tas,'sense',tagDict)
lDict_tas_antisense, tagDict = basicProcessModule_slim.getOverlappingLoci_v02(ghc_all,lc_tas,agi_tas,'antisense',tagDict)

print 'Number of tags:',len(tagDict.keys())


#print contents of each lDict, a modified gHitsAll (with gene and gene-type tags) and separate ghcs

print 'Printing lDicts...'

cmd = 'mkdir -p ' + directory + '/lDicts/'
os.system(cmd)

printLdict(lDict_te_sense,'lDict_te_sense')
printLdict(lDict_te_antisense,'lDict_te_antisense')

printLdict(lDict_mi_sense,'lDict_mi_sense')
printLdict(lDict_mi_antisense,'lDict_mi_antisense')

printLdict(lDict_tas_sense,'lDict_tas_sense')
printLdict(lDict_tas_antisense,'lDict_tas_antisense')


#print lDicts that match spike-ins to outFile

lDict_sRNAspikes_sense,lDict_sRNAspikes_antisense = basicProcessModule_slim.getOverlappingSpikeHits_v02(ghc_sRNAspike)
printLdict_spikes(lDict_sRNAspikes_sense,'lDict_sRNAspikes_sense')
printLdict_spikes(lDict_sRNAspikes_antisense,'lDict_sRNAspikes_antisense')


#print GHCs for TEs, miRNAs and tasiRNAs

print 'Printing ghcs...'

cmd = 'mkdir -p ' + directory + '/gHits_all/'
os.system(cmd)

cmd = 'mkdir -p ' + directory + '/gHits_singleGeneTypes/'
os.system(cmd)

nameList = ['te','mi','tas']

ofh_dict = {}

for name in nameList:
    outFile1 = directory + '/gHits_all/' + name + '_gHits'
    ofh1 = open(outFile1,'w')
    print >>ofh1,'name','\t','seq','\t','chromo','\t','start','\t','end','\t','strand','\t','readNum_raw','\t','readNum_hits','\t','readNum_rpm','\t','geneList','\t','typeList'
    outFile2 = directory + '/gHits_singleGeneTypes/' + name + '_gHits'
    ofh2 = open(outFile2,'w')
    ofh_dict[name] = (ofh1,ofh2)
    print >>ofh2,'name','\t','seq','\t','chromo','\t','start','\t','end','\t','strand','\t','readNum_raw','\t','readNum_hits','\t','readNum_rpm','\t','geneList','\t','typeList'

readTotal_single = 0
readTotal_multi = 0
readTotal_unAnno = 0
readTotal_anno = 0

te_sense_reads = 0
te_antisense_reads = 0
mi_sense_reads = 0
mi_antisense_reads = 0
tas_sense_reads = 0
tas_antisense_reads = 0

te_sense_reads_s = 0
te_antisense_reads_s = 0
mi_sense_reads_s = 0
mi_antisense_reads_s = 0
tas_sense_reads_s = 0
tas_antisense_reads_s = 0

tagList = ghc_all.getHits()

print 'Number of tags in ghc:',len(tagList)

for tag in tagList:
    geneList = tagDict[tag.name()][0]
    typeList = tagDict[tag.name()][1]
    typeList2 = tagDict[tag.name()][2]
    geneList_str = ''
    typeList_str = ''
    for gene in geneList:
        geneList_str = geneList_str + '$' + gene
    for geneType in typeList:
        typeList_str = typeList_str + '$' + geneType
    if len(typeList2) == 0:
        readTotal_unAnno = readTotal_unAnno + tag.readNum_normHit()
    else:
        readTotal_anno = readTotal_anno + tag.readNum_normHit()
        for name in nameList:
            if name in typeList2:
                print_gHit(tag,name)
        for geneType in typeList:
            if geneType == 'te_sense':
                te_sense_reads = te_sense_reads + tag.readNum_normHit()
            elif geneType == 'te_antisense':
                te_antisense_reads = te_antisense_reads + tag.readNum_normHit()  
            elif geneType == 'mi_sense':
                mi_sense_reads = mi_sense_reads + tag.readNum_normHit()
            elif geneType == 'mi_antisense':
                mi_antisense_reads = mi_antisense_reads + tag.readNum_normHit()
            elif geneType == 'tas_sense':
                tas_sense_reads = tas_sense_reads + tag.readNum_normHit()
            elif geneType == 'tas_antisense':
                tas_antisense_reads = tas_antisense_reads + tag.readNum_normHit()
        if len(typeList2) == 1:
            readTotal_single = readTotal_single + tag.readNum_normHit()
            for geneType in typeList: 
                if geneType == 'te_sense':
                    te_sense_reads_s = te_sense_reads_s + tag.readNum_normHit()
                elif geneType == 'te_antisense':
                    te_antisense_reads_s = te_antisense_reads_s + tag.readNum_normHit()  
                elif geneType == 'mi_sense':
                    mi_sense_reads_s = mi_sense_reads_s + tag.readNum_normHit()
                elif geneType == 'mi_antisense':
                    mi_antisense_reads_s = mi_antisense_reads_s + tag.readNum_normHit()
                elif geneType == 'tas_sense':
                    tas_sense_reads_s = tas_sense_reads_s + tag.readNum_normHit()
                elif geneType == 'tas_antisense':
                    tas_antisense_reads_s = tas_antisense_reads_s + tag.readNum_normHit()
        elif len(typeList2) > 1:
            readTotal_multi = readTotal_multi + tag.readNum_normHit()
                

for key in ofh_dict.keys():
    ofh1 = ofh_dict[key][0]
    ofh1.close()
    ofh2 = ofh_dict[key][1]
    ofh2.close()

#print spike gHits

sRNAspike_sense_reads = 0
sRNAspike_antisense_reads = 0

outFile_sS_s = directory + '/gHits_all/sRNAspike_sense_gHits'
outFile_sS_a = directory + '/gHits_all/sRNAspike_antisense_gHits'

ofh_sS_s = open(outFile_sS_s,'w')
ofh_sS_a = open(outFile_sS_a,'w')

for hit in ghc_sRNAspike.getHits():
    if hit.strand() == '+':
        sRNAspike_sense_reads = sRNAspike_sense_reads + hit.readNum_normHit()
        print >>ofh_sS_s,hit.name(),'\t',hit.seq(),'\t',hit.chromo(),'\t',hit.start(),'\t',hit.end(),'\t',hit.strand(),'\t',hit.readNum_raw(),'\t',hit.readNum_normHit(),'\t',hit.readNum_normRpm()
    if hit.strand() == '-':
        sRNAspike_antisense_reads = sRNAspike_antisense_reads + hit.readNum_normHit()
        print >>ofh_sS_a,hit.name(),'\t',hit.seq(),'\t',hit.chromo(),'\t',hit.start(),'\t',hit.end(),'\t',hit.strand(),'\t',hit.readNum_raw(),'\t',hit.readNum_normHit(),'\t',hit.readNum_normRpm()

#check contents

total_reads = ghc_all.get_readNum_normHit()

print 'Total read class','\t','\t','\t','\t','\t','Total number of reads','\t','Percent of total','\t','Percent of annotated'
print 'Full genomic hit collection:','\t','\t','\t','\t',total_reads,'\t','\t','NA','\t','\t','\t','NA'
print 'Reads that do not map to annotated region:','\t','\t',readTotal_unAnno,'\t','\t',round(readTotal_unAnno/float(total_reads),2),'\t','\t','\t','NA'
print 'Reads that map to annotated region:','\t','\t','\t',readTotal_anno,'\t','\t',round(readTotal_anno/float(total_reads),2),'\t','\t','\t','NA'
print 'Total reads that map to multiple anno gene_types:','\t',readTotal_multi,'\t','\t',round(readTotal_multi/float(total_reads),2),'\t','\t','\t',round(readTotal_multi/float(readTotal_anno),2)
print 'Total reads that map to single anno gene_type:','\t','\t',readTotal_single,'\t','\t',round(readTotal_single/float(total_reads),2),'\t','\t','\t',round(readTotal_single/float(readTotal_anno),2)

if sRNA_spikes == 'yes':
    sRNAspike_reads = ghc_sRNAspike.get_readNum_normHit()
    print 'Total reads that map to sRNA spike-ins:','\t','\t',sRNAspike_reads,'\t','\t',round(sRNAspike_reads/float(total_reads + sRNAspike_reads),2)

#now look at individual contents
te_reads = te_sense_reads + te_antisense_reads
mi_reads = mi_sense_reads + mi_antisense_reads
tas_reads = tas_sense_reads + tas_antisense_reads

te_reads_s = te_sense_reads_s + te_antisense_reads_s
mi_reads_s = mi_sense_reads_s + mi_antisense_reads_s
tas_reads_s = tas_sense_reads_s + tas_antisense_reads_s

print '\n'
print 'Gene Type','\t','\t','Total number of reads','\t','Percent of total annotated','\t','Number of reads mapping to single geneType','\t','Percent of single annotated','\t','Percent sense','\t','Percent antisense'
if te_reads_s > 0:
    print 'Transposons:','\t','\t',round(te_reads,2),'\t','\t',round(te_reads/readTotal_anno,2),'\t','\t','\t','\t',round(te_reads_s,2),'\t','\t','\t','\t','\t',round(te_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(te_sense_reads_s/te_reads_s,2),'\t','\t',round(te_antisense_reads_s/te_reads_s,2)
else:
    print 'Transposons:','\t','\t',round(te_reads,2),'\t','\t',round(te_reads/readTotal_anno,2)
if mi_reads_s > 0:
    print 'MIRNA genes:','\t','\t',round(mi_reads,2),'\t','\t',round(mi_reads/readTotal_anno,2),'\t','\t','\t','\t',round(mi_reads_s,2),'\t','\t','\t','\t','\t',round(mi_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(mi_sense_reads_s/mi_reads_s,2),'\t','\t',round(mi_antisense_reads_s/mi_reads_s,2)
else:
    print 'MIRNA genes:','\t','\t',round(mi_reads,2),'\t','\t','\t',round(mi_reads/readTotal_anno,2)
if tas_reads_s > 0:
    print 'TAS genes:','\t','\t',round(tas_reads,2),'\t','\t',round(tas_reads/readTotal_anno,2),'\t','\t','\t','\t',round(tas_reads_s,2),'\t','\t','\t','\t','\t',round(tas_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(tas_sense_reads_s/tas_reads_s,2),'\t','\t',round(tas_antisense_reads_s/tas_reads_s,2)
else:
    print 'TAS genes:','\t','\t',round(tas_reads,2),'\t','\t',round(tas_reads/readTotal_anno,2)

print '#note: total of all gene type categories reads mapping to single annotation gene-type is > total of reads mapping to single anno gene_types because some reads map with same gene-type (e.g. protein-coding genes), but to different strands (e.g. sense vs. antisense)'

#open summary file and print

outFile = directory + '/alignSummary'
ofh = open(outFile,'a')

print >>ofh,'Number of tags in ghc:','\t',len(tagList)
print >>ofh,'Total read class','\t','\t','\t','\t','\t','Total number of reads','\t','Percent of total','\t','Percent of annotated'
print >>ofh,'Full genomic hit collection:','\t','\t','\t','\t',total_reads,'\t','\t','NA','\t','\t','\t','NA'
print >>ofh,'Reads that do not map to annotated region:','\t','\t',readTotal_unAnno,'\t','\t',round(readTotal_unAnno/float(total_reads),2),'\t','\t','\t','NA'
print >>ofh,'Reads that map to annotated region:','\t','\t','\t',readTotal_anno,'\t','\t',round(readTotal_anno/float(total_reads),2),'\t','\t','\t','NA'
print >>ofh,'Total reads that map to multiple anno gene_types:','\t',readTotal_multi,'\t','\t',round(readTotal_multi/float(total_reads),2),'\t','\t','\t',round(readTotal_multi/float(readTotal_anno),2)
print >>ofh,'Total reads that map to single anno gene_type:','\t','\t',readTotal_single,'\t','\t',round(readTotal_single/float(total_reads),2),'\t','\t','\t',round(readTotal_single/float(readTotal_anno),2)

if sRNA_spikes == 'yes':
    sRNAspike_reads = ghc_sRNAspike.get_readNum_normHit()
    print >>ofh,'Total reads that map to sRNA spike-ins:','\t','\t',sRNAspike_reads,'\t','\t',round(sRNAspike_reads/float(total_reads + sRNAspike_reads),2)


#now look at individual contents

print >>ofh,'\n'
print >>ofh,'Gene Type','\t','\t','Total number of reads','\t','Percent of total annotated','\t','Number of reads mapping to single geneType','\t','Percent of single annotated','\t','Percent sense','\t','Percent antisense'
if te_reads_s > 0:
    print >>ofh,'Transposons:','\t','\t',round(te_reads,2),'\t','\t',round(te_reads/readTotal_anno,2),'\t','\t','\t','\t',round(te_reads_s,2),'\t','\t','\t','\t','\t',round(te_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(te_sense_reads_s/te_reads_s,2),'\t','\t',round(te_antisense_reads_s/te_reads_s,2)
else:
    print >>ofh,'Transposons:','\t','\t',round(te_reads,2),'\t','\t',round(te_reads/readTotal_anno,2)
if mi_reads_s > 0:
    print >>ofh,'MIRNA genes:','\t','\t',round(mi_reads,2),'\t','\t',round(mi_reads/readTotal_anno,2),'\t','\t','\t','\t',round(mi_reads_s,2),'\t','\t','\t','\t','\t',round(mi_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(mi_sense_reads_s/mi_reads_s,2),'\t','\t',round(mi_antisense_reads_s/mi_reads_s,2)
else:
    print >>ofh,'MIRNA genes:','\t','\t',round(mi_reads,2),'\t','\t','\t',round(mi_reads/readTotal_anno,2)
if tas_reads_s > 0:
    print >>ofh,'TAS genes:','\t','\t',round(tas_reads,2),'\t','\t',round(tas_reads/readTotal_anno,2),'\t','\t','\t','\t',round(tas_reads_s,2),'\t','\t','\t','\t','\t',round(tas_reads_s/readTotal_single,2),'\t','\t','\t','\t',round(tas_sense_reads_s/tas_reads_s,2),'\t','\t',round(tas_antisense_reads_s/tas_reads_s,2)
else:
    print >>ofh,'TAS genes:','\t','\t',round(tas_reads,2),'\t','\t',round(tas_reads/readTotal_anno,2)

print >>ofh,'#note: total of all gene type categories reads mapping to single annotation gene-type is > total of reads mapping to single anno gene_types because some reads map with same gene-type (e.g. protein-coding genes), but to different strands (e.g. sense vs. antisense)'


ofh.close()

print 'Job complete!'
