#analyzeSmallRNAspikeIns_justDoseResponse.py
#performs basic analysis of small RNA spike-ins
#Input:
#   1. sample
#   2. dilution = dilution factor used to account for amount of RNA used out of total
#Output:
#   1. table of small RNA spike-in name, sRNA-Seq reads per million genome-matching reads (RPM), molar amount and number of molecules

import sys, os, ghcModule_slim

sample = sys.argv[1]
dilution = float(sys.argv[2])       #in this study, all dilutions were 0.0050; 0.0050 = (1 ul of sRNA spike-in mix) * (1/100 dilution of sRNA spike-in prior to adding) * (1/0.5 ug total RNA) * (1/4 ul of sRNA spike-in mix) 

#returns genomic hit collections (one for each spike-in parent) for reads that overlap small RNA spike-ins

def createGHCs_sRNAspike(inFile):
    winSize = 500
    hitList_361 = []
    hitList_403 = []
    hitList_433 = []
    hitList_71 = []
    hitList_823 = []
    hitList_87 = []
    hitList_871 = []
    hitList_974 = []
    hitList_all = []
    ifh = open(inFile)
    inLines = ifh.readlines()
    ifh.close()

    for line in inLines:
        cols = line.strip().split('\t')
        tag = cols[0].strip()
        tagSeq = cols[1].strip()
        chromo = 'cmiR' + cols[2].strip().split('.')[0]
        start = int(cols[3].strip())
        end = int(cols[4].strip())
        strand = cols[5].strip()
        readNum_raw = float(cols[6].strip())
        readNum_normHit = float(cols[7].strip())
        readNum_normRpm = float(cols[8].strip())
        gHit = ghcModule_slim.GenomicHit_rpm(tag,tagSeq,chromo,start,end,strand,readNum_raw,readNum_normHit,readNum_normRpm)
        hitList_all.append(gHit)
        if chromo == 'cmiR361':
            hitList_361.append(gHit)
        if chromo == 'cmiR403':
            hitList_403.append(gHit)
        if chromo == 'cmiR433':
            hitList_433.append(gHit)
        if chromo == 'cmiR71':
            hitList_71.append(gHit)
        if chromo == 'cmiR823':
            hitList_823.append(gHit)
        if chromo == 'cmiR87':
            hitList_87.append(gHit)
        if chromo == 'cmiR871':
            hitList_871.append(gHit)
        if chromo == 'cmiR974':
            hitList_974.append(gHit)

    if len(hitList_all) > 0:
        ghc_all = ghcModule_slim.GenomicHitCollection_rpm(hitList_all,winSize)
    else:
        ghc_all = '.'
    if len(hitList_361) > 0:
        ghc_361 = ghcModule_slim.GenomicHitCollection_rpm(hitList_361,winSize)
    else:
        ghc_361 = '.'
    if len(hitList_403) > 0:
        ghc_403 = ghcModule_slim.GenomicHitCollection_rpm(hitList_403,winSize)
    else:
        ghc_403 = '.'
    if len(hitList_433) > 0:
        ghc_433 = ghcModule_slim.GenomicHitCollection_rpm(hitList_433,winSize)
    else:
        ghc_433 = '.'
    if len(hitList_71) > 0:
        ghc_71 = ghcModule_slim.GenomicHitCollection_rpm(hitList_71,winSize)
    else:
        ghc_71 = '.'
    if len(hitList_823) > 0:
        ghc_823 = ghcModule_slim.GenomicHitCollection_rpm(hitList_823,winSize)
    else:
        ghc_823 = '.'
    if len(hitList_87) > 0:
        ghc_87 = ghcModule_slim.GenomicHitCollection_rpm(hitList_87,winSize)
    else:
        ghc_87 = '.'
    if len(hitList_871) > 0:
        ghc_871 = ghcModule_slim.GenomicHitCollection_rpm(hitList_871,winSize)
    else:
        ghc_871 = '.'
    if len(hitList_974) > 0:
        ghc_974 = ghcModule_slim.GenomicHitCollection_rpm(hitList_974,winSize)
    else:
        ghc_974 = '.'

    return ghc_all,ghc_433,ghc_403,ghc_361,ghc_871,ghc_974,ghc_71,ghc_823,ghc_87

def getCmiRquants(dilution, mixType):
    inFile = 'Ath_annotations/quantTable_' + mixType
    ifh = open(inFile)
    inLines = ifh.readlines()

    quantDict = {}

    for line in inLines:
        if line.startswith('cmi'):
            cols = line.split('\t')
            name = cols[0].strip()
            amol = float(cols[1].strip()) * dilution
            molNum = float(cols[2].strip()) * 1000000 * dilution
            quantDict[name] = (amol,molNum)

    return quantDict


workDir1 = sample + '/lDicts'
workDir2 = sample + '/gHits_all'

#retrieve lDict and gHit files

fileList1 = os.listdir(workDir1)
fileList2 = os.listdir(workDir2)

ldCount = 0
ghCount = 0

for file in fileList1:
    if file.strip() == 'lDict_sRNAspikes_sense':
        ldFile = workDir1 + '/' + file
        ldCount += 1

for file in fileList2:
    if file.strip() == 'sRNAspike_sense_gHits':
        ghFile = workDir2 + '/' + file
        ghCount += 1

if ldCount == 1:
    print 'OK: only one sRNA locusDict file found'
else:
    print 'Caution! More than one sRNA locusDict file found'

if ghCount == 1:
    print 'OK: only one sRNA gHit file found'
else:
    print 'Caution! More than one sRNA gHit file found'


#make neccesary directories

outDir = sample + '/smRNA_spikeIns/'
cmd = 'mkdir -p ' + outDir
os.system(cmd)


#build ghcs for spike-in reads

ghcList = createGHCs_sRNAspike(ghFile)

ghcCount = 0

for ghc in ghcList:
    if ghc != '.':
        ghcCount += 1

print 'Number of ghcs built for cmiRXXX spike-ins:','\t',ghcCount - 1


#calculate rpm and se for each set of cmiRXXX spike-ins, print dose-response curves

outFile = outDir + 'doseResponseTable_noTransform'
ofh = open(outFile,'w')


##retrieve estimated amol/molecules per cmiRXXX values

quantDict = getCmiRquants(dilution, '140126')

print >>ofh,'cmiR','\t','rpm (mean)','\t','amol','\t','mlcs'


#calculate rpm for each parent

for ghc in ghcList[1:]:
    if ghc != '.':
        gHits = ghc.getHits()
        rpm = 0
        for hit in gHits:
            parent = hit.chromo()
            rpm = rpm + hit.readNum_normRpm()
        mlc = quantDict[parent][1]
        print >>ofh,parent,'\t',rpm,'\t',quantDict[parent][0],'\t',quantDict[parent][1]

ofh.close()
