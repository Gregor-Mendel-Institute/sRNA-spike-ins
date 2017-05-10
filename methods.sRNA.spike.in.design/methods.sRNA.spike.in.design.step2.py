#methods.sRNA.spike.in.design.step2.py
#
#input1: noMatch.fa (randomly generated 13mer sequences from step 1 that do not match the genome-of-interest
#           ; Arabidopsis thaliana Col-0 for this study)
#output1: sets of 65,536 possible sequences for each of the above 13mers + 4 random bases on 5'
#           and 3' ends
#
#input2: output1
#output2: 1) folding structures calculated by RNAfold (in folded directory), 2) minimum free energies (MFEs) calculated by RNAfold
            #(in mfes directory) and 3) summary statistics for MFEs (minimum, maximum, mean and median) for each series of random
            #13mers that do not match the genome, and which have 256 possible random 4 bases added to both sides

import random, sys, numpy, os

#functions

def generateBaseList():
    baseList = ['A','C','G','T']
    return baseList

def generateRandomers(inList):
    outList = []
    baseList = generateBaseList()
    while len(baseList) > 0:
        base = baseList.pop()
        for seq in inList:        
            seq = seq + base
            outList.append(seq)
    return outList


inFile = 'noMatch.fa'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

faDict = {}

for i in range(0,len(inLines),2):
    name = inLines[i].strip()
    name = name.lstrip('>')
    seq = inLines[i+1].strip()
    faDict[name] = (seq)

for name in faDict.keys():
    inList = generateBaseList()
    random2mers = generateRandomers(inList)
    random3mers = generateRandomers(random2mers)
    random4mers = generateRandomers(random3mers)
    randomDict = {}
    i = 0
    dict5prime = {}
    seq = faDict[name]
    for mer in random4mers:
        newSeq = mer + seq
        newKey = name + '.' + str(i)
        dict5prime[newKey] = (newSeq)
        i += 1
    i = 0
    for key in dict5prime.keys():
        seq = dict5prime[key]
        for mer in random4mers:
            newSeq = seq + mer
            newKey = key.split('.')
            newKey = newKey[0]
            newKey = newKey + '.' + str(i)
            randomDict[newKey] = (newSeq)
            i += 1

    outFile = 'randomOligoSets/' + name + '.fa'

    ofh = open(outFile,'w')

    for key in randomDict.keys():
        print >>ofh,'>' + key
        print >>ofh,randomDict[key]

    ofh.close()

#fold each set of random oligos and calculate min, max, mean and median
##note: the same commands were used to fold and 

workDir = 'randomOligoSets/'

outDir = 'randomOligoSets/folded/'

fileList = os.listdir(workDir)

for file in fileList:
    if file.endswith('.fa'):
        cmd = 'cat ' + workDir + file + '| RNAfold -T 4 --noPS > ' + outDir + file + '_folded'

        os.system(cmd)

#retrieve mfes, calculate min, max, mean and median, and print to outfiles (summary stats in one and list of mfes for each randomer in the other)

workDir = 'randomOligoSets/folded/'

outDir = 'randomOligoSets/folded/mfes/'

fileList = os.listdir(workDir)

outFile1 = outDir + 'summaryStats.xls'

ofh1 = open(outFile1,'w')

print >>ofh1,'oligo','\t','minimum','\t','maximum','\t','mean','\t','median'

for oligo in faDict.keys():

    outFile2 = outDir + oligo + '_mfes'

    ofh2 = open(outFile2,'w')

    for file in fileList:
        if file.endswith(oligo + '.fa_folded'):
            inFile = workDir + file
            ifh = open(inFile)
            inLines = ifh.readlines()
            ifh.close()

            mfeList = []
            
            for line in inLines:
                if line.startswith('.') or line.startswith('(') or line.startswith(')'):
                    line = line.strip()
                    mfe = float(line[-7:-1])
                    mfeList.append(mfe)

            for mfe in mfeList:
                print >>ofh2,mfe
            
            minimum = round(numpy.min(mfeList),2)
            maximum = round(numpy.max(mfeList),2)
            mean = round(numpy.mean(mfeList),2)
            median = round(numpy.median(mfeList),2)

            print >>ofh1,oligo,'\t',minimum,'\t',maximum,'\t',mean,'\t',median

ofh1.close()
ofh2.close()
    

