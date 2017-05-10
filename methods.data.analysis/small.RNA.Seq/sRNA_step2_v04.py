#sRNA_step2_v02.py
#updated on 17.04.18
#intended to be used after sRNA_step1_var_qual.py
#converts FASTQ file with multiple reads per tag (i.e. unique sequence) to
#FASTA file with single reads per tag, and the number of reads per tag indicated
#beside tag name and separated by an underscore
#Takes folder containing *trimmed.fastq file from sRNA_step1_var_qual.py as input,
#writes tracking information to tracking file contained within the same folder, and
#writes output FA file to same folder

import os

import sys

dataRoot = '/lustre/scratch/users/michael.nodine/seq/'
genomeRoot = '/lustre/scratch/users/michael.nodine/compFiles/genomes/Ath/Col/'


def make1024_dicts():
    parentDict = {}
    baseList = ['A','T','C','G']
    for base1 in baseList:
        string = ''
        for base2 in baseList:
            for base3 in baseList:
                for base4 in baseList:
                    for base5 in baseList:
                        string = base1 + base2 + base3 + base4 + base5
                        parentDict[string] = ({})
    return parentDict


def parse_first5nt_dict(tag,parentDict):
    index = tag[0:5]
    daughterDict = parentDict[index]

    if tag in daughterDict.keys():
        hitCount = daughterDict[tag] + 1
        daughterDict[tag] = (hitCount)
    else:
        hitCount = 1
        daughterDict[tag] = (hitCount)

    parentDict[index] = daughterDict
    
    return parentDict


folder = sys.argv[1]
workDir = dataRoot + folder
fileList = os.listdir(workDir)
trFile = workDir + '/tracking'
trfh = open(trFile,'a')

print >>trfh,'\n'
print >>trfh,'Currently processing trimmed fastq file through sRNA_step2_v02.py'

fileCount = 0

for file in fileList:
    if file.endswith('_trimmed.fastq'):
        fastq = file
        outFileList = file.split('.fastq')
        outFile = workDir + '/' + outFileList[0] + '.fa'
        fileCount += 1

if fileCount > 1:
    print >>trfh,'Error: more than one _trimmed.fastq file detected'
else:
    print >>trfh,'OK: only one _trimmed.fastq file detected'

fastqFile = workDir + '/' + fastq

print >>trfh,'fastqFile:',fastqFile

totalTagNum = 0
a = 0
b = 4


#make all parentDict of 1024 possible subDicts each of which contain seqs with different first four nts

parentDict = make1024_dicts()

print >>trfh,'Currently condensing reads....'


#sort seqs into appropriate daughter dicts keeping track of number of occurences (i.e. reads)

fastqFile = open(fastqFile)
lastLine = False
while lastLine == False:
    package = []
    for i in range(a, b):
        infileLine = fastqFile.readline()
        if infileLine == '':
            lastLine = True
        else:
        ## store line in package
            lastLine = False
            package.append(infileLine.rstrip())
        a = a + 4
        b = b + 4
    if lastLine != True:
        seq = package[1]
        totalTagNum = totalTagNum + 1
        #send to appropriate daughter dict and update
        parentDict = parse_first5nt_dict(seq,parentDict)


#print seqs to outFile

print >>trfh,'Printing .fa file...'

ofh = open(outFile,'w')

i = 0

for key in parentDict.keys():
    dDict = parentDict[key]
    if len(dDict.keys()) > 0:
        for tag in dDict.keys():
            readNum = dDict[tag]
            output = '>' + str(i) + '_' + str(readNum)
            print >>ofh,output
            print >>ofh,tag
            i += 1

ofh.close()

print >>trfh,'.fa file has been printed to ' + outFile
