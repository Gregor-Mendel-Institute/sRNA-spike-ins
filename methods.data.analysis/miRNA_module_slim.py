#miRNA_module.py
#contains functions that deal with miRNA and tasiRNA annotations

import sys, ghcModule_slim

annoRoot = '/lustre/scratch/users/michael.nodine/apps/Ath_annotations/'      #this will have to be set to the appropriate directory

winSize = 500

dataRoot = '/lustre/scratch/users/michael.nodine/apps/Ath_annotations/seq/sRNA/'

#15.02.17: used to search for ghc that overlap MIRNA

def getOverlappingLoci_miRNA(ghc,lc,strand):
    hitDict = {} #[miRNA] = (gHitList)
    lc_list = lc.getLoci()

    for locus in lc_list: 
        overlaplist = ghc.getOverlap(locus,strand)  #get gHits that overlap locus
        for j in overlaplist:                       #for each gHit that overlaps locus:
            if locus in hitDict.keys():             #if already in hitDict_genomic keys, then append to locusList      
                locusList = hitDict[locus]
                locusList.append(j)
                hitDict[locus] = (locusList)
            else:                                       #else add as new entry
                locusList = [j]
                hitDict[locus] = (locusList)
            
    return hitDict


#takes matList in format given in mir class and returns as list of loci

def getMatList(mir):
    locusList = []
    
    matList = mir.matList()

    m = matList
    #print m
    mList = str(m).split(',')
    for m in mList:
        m = m.strip()
        m = m.lstrip("['")
        m = m.strip("\n']")
        m = m.lstrip(" '")
        m = m.split('$')
        name = m[0]
        chromo = m[1]
        start = m[2]
        end = m[3]
        strand = m[4][0]
        fam = m[5]
        locus = ghcModule_slim.Locus(name,chromo,start,end,strand,'miRNA',mir.common(),fam,'')
        locusList.append(locus)

    return locusList

def getMatList_tasi(tas):
    locusList = []
    
    matList = tas.long_desc()
    #print tas.short_desc()
    #print 'matList:',matList
    
    #if len(matList) > 0:
        #print len(matList)
    for mat in matList:
        subList = mat.split('$')
        #print 'subList:',subList
        name = subList[0]
        chromo = subList[1]
        start = int(subList[2])
        end = int(subList[3])
        strand = subList[4]

        locus = ghcModule_slim.Locus(name,chromo,start,end,strand,'tasiRNA',tas.short_desc(),'','')
        locusList.append(locus)

    return locusList


#returns ghc from list gHits (e.g. hits that overlap MIRNA genes)

def makeGHCfromHitList(hitList):
    gList = []
    
    for hit in hitList:
        gHit = ghcModule_slim.GenomicHit_rpm(hit.name(),hit.seq(),hit.chromo(),hit.start(),hit.end(),hit.strand(),hit.readNum_raw(),hit.readNum_normHit(),hit.readNum_normRpm())
        gList.append(gHit)

    ghc = ghcModule_slim.GenomicHitCollection_rpm(gList,winSize)

    return ghc


#returns gHits (in form of ghc) that are contained within given loci with given extension

def getMirContain(locus,ghc,extension):
    locus = ghcModule_slim.Locus(locus.geneID(),locus.chromo(),locus.start() - extension,locus.end() + extension,locus.strand(),'','','','')
    containList = []
    for gHit in ghc.getHits():
        hitLocus = ghcModule_slim.Locus(gHit.name(),gHit.chromo(),gHit.start(),gHit.end(),gHit.strand(),'','','','')
        if locus.contains(hitLocus) == True:
            containList.append(gHit)

    ghc = ghcModule_slim.GenomicHitCollection(containList,winSize)

    return ghc


def getTasContain(locus,ghc,extension):
    locus = ghcModule_slim.Locus(locus.geneID(),locus.chromo(),locus.start() - extension,locus.end() + extension,locus.strand(),'','','','')
    containList = []
    for gHit in ghc.getHits():
        hitLocus = ghcModule_slim.Locus(gHit.name(),gHit.chromo(),gHit.start(),gHit.end(),gHit.strand(),'','','','')
        #print hitLocus.chromo(),hitLocus.start(),hitLocus.end(),hitLocus.strand()
        #print locus.contains(hitLocus)
        if locus.contains(hitLocus) == True:
            gHit_tas = ghcModule_slim.GenomicHit_rpm(locus.geneID(),gHit.seq(),gHit.chromo(),gHit.start(),gHit.end(),gHit.strand(),gHit.readNum_raw(),gHit.readNum_normHit(),gHit.readNum_normRpm())
            containList.append(gHit_tas)

    ghc = ghcModule_slim.GenomicHitCollection(containList,winSize)

    return ghc


#looks for annotated mature sequences and returns overlapping tags (+/- 2 nt)

def getContained_mature(hitDict):
    extension = 2
    matDict = {} #[mature_locus] = (ghc) ; [mature miRNA as locus class] = (genomic hit collection of all reads contained within)
    
    for mir in hitDict.keys():
        locusList = getMatList(mir)
        hitList = hitDict[mir]
        ghc = makeGHCfromHitList(hitList)
        for locus in locusList:
            ghc_mir = getMirContain(locus,ghc,extension)
            matDict[locus] = (ghc_mir)

    return matDict

def getContained_tasi(hitDict):
    extension = 2
    matDict = {} #[mature_locus] = (ghc) ; [mature tasiRNA as locus class] = (genomic hit collection of all reads contained within)
    
    for tasi in hitDict.keys():
        locusList = getMatList_tasi(tasi)
        hitList = hitDict[tasi]
        ghc = makeGHCfromHitList(hitList)
        for locus in locusList:
            ghc_mir = getTasContain(locus,ghc,extension)
            if len(ghc_mir.getHits()) > 0:
                matDict[locus] = (ghc_mir)

    return matDict

#returns the most abundant sequence from a given ghc

def getMostAbundFromGHC(ghc):
    maxRpm = 0
    maxHit = ''
    gHits = ghc.getHits()
    for hit in gHits:
        if hit.readNum_raw() > maxRpm:
            maxRpm = hit.readNum_raw()
            maxHit = hit

    return maxHit


#returns most abundant gHit for all reads overlapping MIRNA genes

def getAbund(hitDict):
    abundDict = {}

    for mir in hitDict.keys():
        hitList = hitDict[mir]
        ghc = makeGHCfromHitList(hitList)
        maxHit = getMostAbundFromGHC(ghc)
        abundDict[mir] = (maxHit)
    return abundDict

def getAnnoFamDict():
    inFile = annoRoot + '/miRNA/annoFams_150828'
    ifh = open(inFile)
    inLines = ifh.readlines()

    famDict = {}

    famList = []
    
    for line in inLines[2:]:
        cols = line.split('\t')
        mir = cols[0].strip()
        fam = cols[1].strip()
        if fam not in famList:
            famList.append(fam)
        if mir in famDict.keys():
            print "ERROR! miRNA entry already exists!"
        else:
            famDict[mir] = (fam)

    print 'Number of annotated miRNAs:',len(famDict.keys())
    print 'Number of annotated miRNA families:',len(famList)

    return famDict

#retrieves dict[tasiRNA] = (fam) where tasiRNAs have Allen et al. nomenclature and fam is if they map
#to tasiRNAdb annotations with three or less mismatches; in this version four TAS1-derived smRNAs were
#grouped into tasiR-HTT family because they are predicted to target Heat-induced TAS1 Targets (HTT) RNAs

def getAnnoFamDict_tasi():
    inFile = annoRoot + '/tasiRNA/annoFams_alt_150921'
    ifh = open(inFile)
    inLines = ifh.readlines()

    famDict = {}
    famList = []

    for line in inLines[2:]:
        cols = line.split('\t')
        mir = cols[0].strip()
        fam = cols[1].strip()
        if fam != '' and fam not in famList:
            famList.append(fam)
        if fam == '':
            fam = mir
        if mir in famDict.keys():
            print "ERROR! tasiRNA entry already exists!"
        else:
            famDict[mir] = (fam)

    print 'Number of annotated tasiRNAs:',len(famDict.keys())
    print 'Number of annotated tasiRNA families:',len(famList)

    return famDict
