#bwtToAllHits_v03.py
#converts genome.bwt file to tab-delimited file (allgHits) with following features:
#name, seq, chromo, start, end, strand, readNum_raw, readNum_hits, readNum_rpm
#also prints summary of alignments and normalizations to alignSummary file
#
#input:
#   1: filePath; path where files associated with given sample are located
#   2: sample; name of sample 
#   3: total_reads; number of genome-matching reads to use for normalization
#

import sys

filePath = sys.argv[1].strip()
sample = sys.argv[2].strip()
total_reads = float(sys.argv[3].strip())/1000000

directory = filePath + '/' + sample + '/'

def makeRC(seq):
    i = 0
    newSeq = ''
    for i in range(0,len(seq)):
        if seq[i] == 'A':
            newSeq = 'T' + newSeq
        if seq[i] == 'T':
            newSeq = 'A' + newSeq
        if seq[i] == 'G':
            newSeq = 'C' + newSeq
        if seq[i] == 'C':
            newSeq = 'G' + newSeq
        if seq[i] == 'N':
            newSeq = 'N' + newSeq
    return newSeq


print "Total number of (million) genome-matching reads:",total_reads

#first make dictionary of all possible tags

inFile = directory + '/tags.genome.bwt'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

hitNumDict = {}

for line in inLines:
    cols = line.split('\t')
    name = cols[0].strip()
    chromo = cols[2].strip()
    hitNumDict[name] = (0)

print "Total number of unique tags:",len(hitNumDict.keys())

#now tally number of hits for each tag

inFile = directory + '/tags.genome.bwt'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

for line in inLines:
    cols = line.split('\t')
    name = cols[0].strip()
    chromo = cols[2].strip()
    hitNum = hitNumDict[name] + 1
    hitNumDict[name] = (hitNum)

print "Length of hitNumDict:",len(hitNumDict.keys())

 
#normalize and print to allgHits outFile

outFile = directory + '/allgHits'
ofh = open(outFile,'w')

print >>ofh,'name','\t','seq','\t','chromo','\t','start','\t','end','\t','strand','\t','readNum_raw','\t','readNum_hits','\t','readNum_rpm'

inFile = directory + '/tags.genome.bwt'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

readNum = 0
hitCount = 0

for line in inLines:
    cols = line.split('\t')
    name = cols[0].strip()
    chromo = cols[2].strip()
    strand = cols[1].strip()
    if strand == '+':
        seq = cols[4].strip()
    elif strand == '-':
        seq = makeRC(cols[4].strip())
    start = int(cols[3].strip())
    end = start + len(seq)
    readNum_raw = int(name.split('_')[1])
    readNum_hits = float(readNum_raw)/(float(hitNumDict[name]))
    readNum_rpm = readNum_hits/total_reads
    readNum = readNum + readNum_hits
    hitCount += 1
    print >>ofh,name,'\t',seq,'\t',chromo,'\t',start,'\t',end,'\t',strand,'\t',readNum_raw,'\t',readNum_hits,'\t',readNum_rpm            

ofh.close()

print 'Number of aligned tags:',len(hitNumDict.keys())
print 'Number alignments:',hitCount
print 'Total (hit-normalized) reads aligned to genome:',readNum,'\t',total_reads

outFile = directory + '/alignSummary'
ofh = open(outFile,'w')

print >>ofh,'\n'
print >>ofh,'Number of aligned tags:','\t',len(hitNumDict.keys())
print >>ofh,'Number alignments:','\t',hitCount
print >>ofh,'Total (hit-normalized) reads aligned to genome (and spike-ins):','\t',readNum,'\t',total_reads

ofh.close()
