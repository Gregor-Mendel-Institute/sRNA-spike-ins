#getLengthVals.py
#generates table of RPM levels for specified range of read lengths; all tallies the RPM for each
#of the four bases at a specified position within the read (e.g. 0 indicates first position)
#input:
#   1. comma-separated list of samples to generate tables for; in this study it was 'col0_leaf1,col0_leaf2,col0_fb1,col0_fb2,d234_fb1,d234_fb2'
#   2. smallest sRNA-Seq read size to be considered; in this study it was '18'
#   3. largest sRNA-Seq read size (+1) to be considered; in this study it was '31'
#   4. type of gene to be considered; in this study it was 'all' meaning to consider all genome-matching reads
#   5. position within read to report at base frequencies; in this study it was '0' indicating to inspect this for the first position


import sys, os

sampleList = sys.argv[1].strip().split(',')
startRange = int(sys.argv[2].strip())
endRange = int(sys.argv[3].strip())
geneType = sys.argv[4].strip()
pos = int(sys.argv[5].strip())

print 'folderList:',sampleList
print 'start of range:',str(startRange)
print 'end of range:',str(endRange)
print 'geneType:',geneType
print 'position for which base frequency is calculated:',str(pos)
    
def getTags(inFile):
    print 'retrieving tags...'
    ifh = open(inFile)
    inLines = ifh.readlines()
    ifh.close()

    tagDict = {}

    for line in inLines[1:]:
        cols = line.strip().split('\t')
        seq = cols[1].strip()
        tag = cols[0].strip()
        rpm = float(cols[8].strip())
        tagDict[tag] = (seq,rpm)

    print 'Number of tags mapping to',geneType,':',len(tagDict.keys()) 
    return tagDict

def retrieveHitDict(directory, geneType):
    if geneType == 'all':
        hitDict_file = directory + '/allgHits'
    else:
        hitDict_file = directory + '/gHits_all/' + geneType + '_gHits'
    
    tagDict = getTags(hitDict_file)      #tagDict[tag]=(length,rpm)

    return tagDict

def getSizeRpm(tagDict):
    print 'retrieving percent of bases at each position...'
    pDict = {}
    for i in range(startRange,endRange):
        pDict[i] = (0,0,0,0)  #pDict[size] = (rpm_A,rpm_G,rpm_C,rpm_T)

    for tag in tagDict.keys():
        seq = tagDict[tag][0]
        size = len(seq)
        if size in range(startRange,endRange):
            #tally = pDict[size][0] + tagDict[tag][1]
            if seq[pos] == 'A':
                tallyA = pDict[size][0] + tagDict[tag][1]
                tallyG = pDict[size][1]
                tallyC = pDict[size][2]
                tallyT = pDict[size][3]
            elif seq[pos] == 'G':
                tallyA = pDict[size][0]
                tallyG = pDict[size][1] + tagDict[tag][1]
                tallyC = pDict[size][2]
                tallyT = pDict[size][3]
            elif seq[pos] == 'C':
                tallyA = pDict[size][0]
                tallyG = pDict[size][1]
                tallyC = pDict[size][2] + tagDict[tag][1]
                tallyT = pDict[size][3]
            elif seq[pos] == 'T':
                tallyA = pDict[size][0]
                tallyG = pDict[size][1]
                tallyC = pDict[size][2]
                tallyT = pDict[size][3] + tagDict[tag][1]
            pDict[size] = (tallyA,tallyG,tallyC,tallyT)

    return pDict


for sample in sampleList:
    print '######################' + sample + '######################'

    outDir = sample + '/lenDistrib/'

    cmd = 'mkdir -p ' + outDir
    print cmd
    os.system(cmd)

    tagDict = retrieveHitDict(sample, geneType)
    pDict = getSizeRpm(tagDict)

    outFile = outDir + 'lenDistrib__geneType_' + geneType + '__start_' + str(startRange) + '__end_' + str(endRange) + '__pos_' + str(pos)
    ofh = open(outFile,'w')
    print >>ofh,'length','\t','A','\t','U','\t','C','\t','G'
    for size in pDict.keys():
        print >>ofh,size,'\t',pDict[size][0],'\t',pDict[size][3],'\t',pDict[size][2],'\t',pDict[size][1]
    ofh.close()

