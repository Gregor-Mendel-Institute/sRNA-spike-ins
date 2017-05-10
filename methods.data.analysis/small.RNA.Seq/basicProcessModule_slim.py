#basicProcessModule_slim.py
#Contains functions required for make_ghcsAndlDicts.py, miR_finder_v03.py and tasiR_finder_anno_v03.py

import os, sys, ghcModule_slim

annoRoot = '/lustre/scratch/users/michael.nodine/apps/Ath_annotations/'      #this will have to be set to the appropriate directory


#makes genomic hit collection (ghc) from allgHits file

def makeGHC_fromAllgHits(filePath,mRNA_spikes,sRNA_spikes):
    winSize = 500
    gHitList = []
    gHitList_mRNA_spikes = []
    gHitList_sRNA_spikes = []

    hitFile = filePath + '/allgHits'
    hfh = open(hitFile)
    inLines = hfh.readlines()
    hfh.close()

    if (mRNA_spikes == 'no') and (sRNA_spikes == 'no'):
        for line in inLines[1:]:
            cols = line.split('\t')
            tag = cols[0].strip()
            seq = cols[1].strip()
            chromo = cols[2].strip()
            start = int(cols[3].strip())
            end = int(cols[4].strip())
            strand = cols[5].strip()
            readNum_raw = int(cols[6].strip())
            readNum_norm = float(cols[7].strip())
            readNum_rpm = float(cols[8].strip())
            gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,strand,readNum_raw,readNum_norm,readNum_rpm)
            gHitList.append(gHit)

        ghc_genomic = ghcModule_slim.GenomicHitCollection_rpm(gHitList,winSize)
        return ghc_genomic

    else:
        for line in inLines[1:]:
            cols = line.split('\t')
            tag = cols[0].strip()
            seq = cols[1].strip()
            chromo = cols[2].strip()
            start = int(cols[3].strip())
            end = int(cols[4].strip())
            strand = cols[5].strip()
            readNum_raw = int(cols[6].strip())
            readNum_norm = float(cols[7].strip())
            readNum_rpm = float(cols[8].strip())

            if chromo.startswith('Ath'):
                gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,strand,readNum_raw,readNum_norm,readNum_rpm)
                gHitList.append(gHit)
            elif mRNA_spikes == 'yes' and chromo.startswith('ERCC'):
                gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,strand,readNum_raw,readNum_norm,readNum_rpm)
                gHitList_mRNA_spikes.append(gHit)
            elif sRNA_spikes == 'yes': # and chromo.startswith('cmiR'):
                gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,strand,readNum_raw,readNum_norm,readNum_rpm)
                gHitList_sRNA_spikes.append(gHit)

        if (mRNA_spikes == 'yes') and (sRNA_spikes == 'yes'):
            ghc_genomic = ghcModule_slim.GenomicHitCollection_rpm(gHitList,winSize)
            ghc_mRNAspike = ghcModule_slim.GenomicHitCollection_rpm(gHitList_mRNA_spikes,winSize)
            ghc_sRNAspike = ghcModule_slim.GenomicHitCollection_rpm(gHitList_sRNA_spikes,winSize)
            return ghc_genomic, ghc_mRNAspike, ghc_sRNAspike

        elif (mRNA_spikes == 'yes'):
            ghc_genomic = ghcModule_slim.GenomicHitCollection_rpm(gHitList,winSize)
            ghc_mRNAspike = ghcModule_slim.GenomicHitCollection_rpm(gHitList_mRNA_spikes,winSize)
            return ghc_genomic, ghc_mRNAspike

        elif (sRNA_spikes == 'yes'):
            ghc_genomic = ghcModule_slim.GenomicHitCollection_rpm(gHitList,winSize)
            ghc_sRNAspike = ghcModule_slim.GenomicHitCollection_rpm(gHitList_sRNA_spikes,winSize)
            return ghc_genomic, ghc_sRNAspike

#makes collection of loci

def makeLocusCollection(locusType):
    #winSize is the size of the denominator used to group loci
    winSize = 500
    if locusType == 'mi':
        annoFile = annoRoot + 'miRBase21_and_TAIR10_miRNA'
    else:
        annoFile = annoRoot + 'TAIR10_' + locusType

    afh = open(annoFile)
    inLines = afh.readlines()

    locusList = []
    lc = None

    if locusType == 'mi':
        for line in inLines[2:]:
            cols = line.split('\t')
            geneID = cols[0].strip()
            chromo = cols[1].strip()
            start = int(cols[2].strip())
            end = int(cols[3].strip())
            strand = cols[4].strip()
            gene_type = locusType
            common = cols[5].strip()
            AGI_ID = cols[6].strip()
            if common == '.':
                common = AGI_ID
            MB_ID = cols[7].strip()
            family = cols[8].strip()
            if family == 'miR156' or family == 'miR157':
                family = 'miR156/157'
            elif family == 'miR165' or family == 'miR166':
                family = 'miR165/166'
            elif family == 'miR170' or family == 'miR171':
                family = 'miR170/171'
            short_desc = cols[9].strip()
            matList = cols[10].strip().split(',')
            pTargets = '.'
            vTargets = '.'
            locus = ghcModule_slim.Locus_miRNA(geneID,chromo,start,end,strand,gene_type,common,AGI_ID,MB_ID,family,short_desc,matList,pTargets,vTargets)
            locusList.append(locus)

        lc = ghcModule_slim.LocusCollection(locusList,winSize)

    elif locusType == 'tas':
        for line in inLines[2:]:
            cols = line.split('\t')
            geneID = cols[0].strip()
            chromo = cols[1].strip()
            start = int(cols[2].strip())
            end = int(cols[3].strip())
            strand = cols[4].strip()
            gene_type = locusType 
            common = cols[6].strip().upper()
            matList = cols[8].strip().split(',')
            long_desc = '.'
            comp_desc = '.'
            locus = ghcModule_slim.Locus(geneID,chromo,start,end,strand,gene_type,common,matList,comp_desc)
            locusList.append(locus)
        lc = ghcModule_slim.LocusCollection(locusList,winSize)
    else:
        for line in inLines[2:]:
            cols = line.split('\t')
            geneID = cols[0].strip()
            chromo = cols[1].strip()
            start = int(cols[2].strip())
            end = int(cols[3].strip())
            strand = cols[4].strip()
            gene_type = locusType 
            if len(cols) > 6:
                short_desc = cols[6].strip()
            else:
                short_desc = '.'
            if len(cols) > 7:
                long_desc = cols[7].strip()
            else:
                long_desc = '.'
            if len(cols) > 8:
                comp_desc = cols[8].strip()
            else:
                comp_desc = '.'
            locus = ghcModule_slim.Locus(geneID,chromo,start,end,strand,gene_type,short_desc,long_desc,comp_desc)
            locusList.append(locus)

        lc = ghcModule_slim.LocusCollection(locusList,winSize)

    return lc


#makeAGIlist makes a list of AGI ids for given type of annotated loci

def makeAGIlist(locusType):
    annoFile = annoRoot + locusType + '_AGIs'
    afh = open(annoFile)
    inLines = afh.readlines()

    agiList = []

    for line in inLines[2:]:
        cols = line.split('\t')
        agi = cols[0].strip()
        agiList.append(agi)
    return agiList


#used to search for gHits that overlap given loci

def getOverlappingLoci_v02(ghc,lc,agiList,strand,tagDict):
    hitDict = {}
    lc_list = lc.getLoci()
    lDict = {}

    for locus in lc_list:
        tagNum = 0
        raw = 0
        hit = 0
        rpm = 0
        tagList = []
        geneType = locus.gene_type() + '_' + strand
        AGI = locus.geneID() + '_' + strand
        overlaplist = []
        overlaplist = ghc.getOverlap(locus,strand)  #get gHits that overlap locus
        if len(overlaplist) > 0:
            for j in overlaplist:                   #for each gHit that overlaps locus:
                tag = j.name()
                tagList.append(tag)
                tagNum = tagNum + 1
                raw = raw + j.readNum_raw()
                hit = hit + j.readNum_normHit()
                rpm = rpm + j.readNum_normRpm()

                geneList = tagDict[tag][0]
                typeList = tagDict[tag][1]
                typeList2 = tagDict[tag][2]
                if AGI not in geneList:
                    geneList.append(AGI)
                if geneType not in typeList:
                    typeList.append(geneType)
                if locus.gene_type() not in typeList2:
                    typeList2.append(locus.gene_type())
                tagDict[tag] = (geneList,typeList,typeList2)

        lDict[locus] = (tagNum,raw,hit,rpm,tagList)

    return lDict, tagDict

def getOverlappingSpikeHits_v02(ghc_spike):
    lDict_sense = {}
    lDict_antisense = {}
    hitList = ghc_spike.getHits()

    for gHit in hitList:
        tagNum = 0
        raw = 0
        hit = 0
        rpm = 0
        tagList = []
        spikeId = gHit.chromo()
        strand = gHit.strand()
        if strand == '+':
            if spikeId in lDict_sense.keys():
                tag = gHit.name()
                tagNum = lDict_sense[spikeId][0] + 1
                readNum_raw = lDict_sense[spikeId][1] + gHit.readNum_raw()
                readNum_hit = lDict_sense[spikeId][2] + gHit.readNum_normHit()
                readNum_rpm = lDict_sense[spikeId][3] + gHit.readNum_normRpm()
                tagList = lDict_sense[spikeId][4]
                tagList.append(tag)
                lDict_sense[spikeId] = (tagNum,readNum_raw,readNum_hit,readNum_rpm,tagList)
            else:
                tag = gHit.name()
                tagList.append(tag)
                tagNum = 1
                lDict_sense[spikeId] = (tagNum,gHit.readNum_raw(),gHit.readNum_normHit(),gHit.readNum_normRpm(),tagList)

        if strand == '-':
            if spikeId in lDict_antisense.keys():
                tag = gHit.name()
                tagNum = lDict_antisense[spikeId][0] + 1
                readNum_raw = lDict_antisense[spikeId][1] + gHit.readNum_raw()
                readNum_hit = lDict_antisense[spikeId][2] + gHit.readNum_normHit()
                readNum_rpm = lDict_antisense[spikeId][3] + gHit.readNum_normRpm()
                tagList = lDict_antisense[spikeId][4]
                tagList.append(tag)
                lDict_antisense[spikeId] = (tagNum,readNum_raw,readNum_hit,readNum_rpm,tagList)
            else:
                tag = gHit.name()
                tagList.append(tag)
                tagNum = 1
                lDict_antisense[spikeId] = (tagNum,gHit.readNum_raw(),gHit.readNum_normHit(),gHit.readNum_normRpm(),tagList)

    return lDict_sense,lDict_antisense


def makeGHC_fromHitsFile_v05(filePath, locusType, collection, strand):
    winSize = 500
    fileList = []
    gHitList = []

    if locusType == 'all':
        inFile = filePath + '/allgHits'
        ifh = open(inFile)
        inLines = ifh.readlines()
        ifh.close()
        for line in inLines[1:]:
            cols = line.split('\t')
            tag = cols[0].strip()
            seq = cols[1].strip()
            chromo = cols[2].strip()
            start = int(cols[3].strip())
            end = int(cols[4].strip())
            tag_strand = cols[5].strip()
            readNum_raw = int(cols[6].strip())
            readNum_norm = float(cols[7].strip())
            readNum_rpm = float(cols[8].strip())
            gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,tag_strand,readNum_raw,readNum_norm,readNum_rpm)
            gHitList.append(gHit)
    else:
        if collection == 'all':
            filePath = filePath + '/gHits_all/'
        else:
            filePath = filePath + '/gHits_singleGeneTypes/'

        if (locusType == 'all_filt'):
            lcsList = ['mi','or','pc','ps','tas','te','unAnno','linc','sORF']
            for lcs in lcsList:
                fileName = lcs + '_gHits'
                fileList.append(fileName)
        else:
            fileName = locusType + '_gHits'
            fileList.append(fileName)

        for fileName in fileList:
            inFile = filePath + '/' + fileName
            ifh = open(inFile)
            inLines = ifh.readlines()
            ifh.close()
            for line in inLines[1:]:
                cols = line.split('\t')
                tag = cols[0].strip()
                seq = cols[1].strip()
                chromo = cols[2].strip()
                start = int(cols[3].strip())
                end = int(cols[4].strip())
                tag_strand = cols[5].strip()
                typeList = cols[10].strip().split('$')
                for lType in typeList:
                    if strand == 'sense' or strand == 'antisense':
                        if lType == locusType + '_' + strand:
                            readNum_raw = int(cols[6].strip())
                            readNum_norm = float(cols[7].strip())
                            readNum_rpm = float(cols[8].strip())
                            gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,tag_strand,readNum_raw,readNum_norm,readNum_rpm)
                            gHitList.append(gHit)
                    elif strand == 'both':
                        if lType == locusType + '_sense' or lType == locusType + '_antisense':
                            readNum_raw = int(cols[6].strip())
                            readNum_norm = float(cols[7].strip())
                            readNum_rpm = float(cols[8].strip())
                            gHit = ghcModule_slim.GenomicHit_rpm(tag,seq,chromo,start,end,tag_strand,readNum_raw,readNum_norm,readNum_rpm)
                            gHitList.append(gHit)
    ghc = ghcModule_slim.GenomicHitCollection_rpm(gHitList,winSize)

    return ghc

