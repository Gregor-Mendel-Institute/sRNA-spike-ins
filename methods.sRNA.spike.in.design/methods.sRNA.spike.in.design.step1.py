#methods.sRNA.spike.in.design.step1.py
#
#input1: mature.miRNA.seqs.top50percent file that contains sequences for most abundant miRNAs
#output1: mature.miRNA.seqs.top50percent.fa
#
#input2: mature.miRNA.seqs.top50percent
#output2: random.fa 

import random

#reads in mature.miRNA.seqs.top50percent file and returns the mature.miRNA.seqs.top50percent.fa file
inFile = 'mature.miRNA.seqs.top50percent'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

outFile = 'mature.miRNA.seqs.top50percent.fa'
ofh = open(outFile,'w')

for line in inLines[2:]:
    if len(line) > 0:
        cols = line.split('\t')
        seq = cols[0].strip()
        name = cols[1].strip()        
        print >>ofh,'>' + name
        print >>ofh,seq

ofh.close()


#reads in mature.miRNA.seqs.top50percent file, and 1) records frequency of each base at each position, 2) based on
#these frequencies generates 1000 random 13mer sub-sequences (positions 4-17)

inFile = 'mature.miRNA.seqs.top50percent'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

freqDict = {}

for i in range(0,21):
    freqDict[i] = (0,0,0,0)

for line in inLines[2:]:
    cols = line.split('\t')
    seq = cols[0].strip()
    for i in range(0,21):
        if len(seq) > i:
            if seq[i] == 'A':
                freqDict[i] = (freqDict[i][0] + 1,freqDict[i][1],freqDict[i][2],freqDict[i][3])
            if seq[i] == 'G':
                freqDict[i] = (freqDict[i][0],freqDict[i][1] + 1,freqDict[i][2],freqDict[i][3])
            if seq[i] == 'C':
                freqDict[i] = (freqDict[i][0],freqDict[i][1],freqDict[i][2] + 1,freqDict[i][3])
            if seq[i] == 'T':
                freqDict[i] = (freqDict[i][0],freqDict[i][1],freqDict[i][2],freqDict[i][3] + 1)

freqDict2 = {}

for pos in freqDict.keys():
    baseList = []
    A_count = freqDict[pos][0]
    G_count = freqDict[pos][1]
    C_count = freqDict[pos][2]
    T_count = freqDict[pos][3]

    for i in range(0,A_count):
        baseList.append('A')
    for i in range(0,G_count):
        baseList.append('G')
    for i in range(0,C_count):
        baseList.append('C')
    for i in range(0,T_count):
        baseList.append('T')
    
    freqDict2[pos] = (baseList)


outFile = 'random.fa'
ofh = open(outFile,'w')

for i in range(0,1000):
    randomSeq = ''
    for pos in freqDict2.keys():
        randomSeq = randomSeq + random.choice(freqDict2[pos])
    print >>ofh,'>' + str(i)
    print >>ofh,randomSeq[4:17]

    
    
        
    
