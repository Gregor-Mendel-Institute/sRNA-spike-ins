#getMirAndTasPrecursorVals.py
#extracts levels (in RPM) for annotated mature miRNAs and tasiRNAs for each corresponding precursor

outDir = "data.for.graphs/"

sampleList = ['col0_leaf1','col0_leaf2','col0_fb1','col0_fb2','d234_fb1','d234_fb2']

def getSrnaVals(inFile,precDict):
    ifh = open(inFile)
    inLines = ifh.readlines()
    ifh.close()

    matDict = {}

    for line in inLines[1:]:
        cols = line.split('\t')
        mat = cols[0].strip()
        rpm = float(cols[1].strip())
        matDict[mat] = (rpm)

    outDict = {}        #[AGI] = (rpm sum of all small RNAs from that precursor)

    for AGI in precDict.keys():
        total_rpm = 0
        sRNAs = precDict[AGI]
        for sRNA in sRNAs:
            if sRNA in matDict.keys():
                total_rpm = total_rpm + matDict[sRNA]
        outDict[AGI] = (total_rpm)

    return outDict

def printOut(outFile, outDict):
    ofh = open(outFile,'w')

    print >>ofh,'AGI','\t','sRNA_rpm'   #,'\t','fpkm_rep2','\t','fpkm_rep3'

    for agi in outDict.keys():
        print >>ofh,agi,'\t',outDict[agi]

    ofh.close()   


#generate dictionaries of precursors and corresponding list of mature sRNAs

precDict_mir = {}           #[AGI] = (miRlist)

inFile = 'Ath_annotations/miRBase21_and_TAIR10_miRNA'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

for line in inLines[2:]:
    cols = line.strip().split('\t')
    embedList = cols[10].strip().split(',')
    miRlist = []
    for item in embedList:
        subList = item.split('$')
        miRlist.append(subList[0])
    AGI = cols[6].strip()
    if AGI == '.':
        AGI = cols[5].strip().upper()

    precDict_mir[AGI] = (miRlist)

print "Number of annotated miRNA precursors:",len(precDict_mir.keys())


precDict_tas = {}
tasAGIlist = []

inFile = 'Ath_annotations/TAIR10_tas'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

for line in inLines[2:]:
    cols = line.strip().split('\t')
    tasList_full = cols[8].strip().split(',')
    tasList = []
    for entry in tasList_full:
        tas = entry.split('$')[0]
        if tas not in tasList:
            tasList.append(tas)
    AGI = cols[0].strip()
    if AGI not in tasAGIlist:
        tasAGIlist.append(AGI)

    precDict_tas[AGI] = (tasList)
            
print "Number of annotated tasiRNA precursors:",len(precDict_tas.keys())
print "Number of TAS RNAs:",len(tasAGIlist)


for sample in sampleList:
    print sample    
    miRfile = sample + '/miRNA/mature_byRpm'
    tasiRfile = sample + '/tasiRNA/mature_byRpm'

    miR_precDict = getSrnaVals(miRfile,precDict_mir)
    tas_precDict = getSrnaVals(tasiRfile,precDict_tas)

    printOut(outDir + sample + '_miR_prec', miR_precDict)
    printOut(outDir + sample + '_tas_prec', tas_precDict)
