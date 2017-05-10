#sRNA_step1_var_qual.py
#updated on 17.04.07
#filters NGS reads in a fastq file based on length of reads, absence of 'N' base calls,
#presence of adaptor sequences, and (optionally) quality of base calls
#input arguments:
##folder: directory in which FASTQ file for NGS reads is located
##lowerLimit: smallest number of bases per NGS reads required; '18' was used in this study
##upperLimit: largest number of bases per NGS reads required; '75' was used in this study
##adapterSeq: adaptor sequence that is required at the 3' end of NGS reads; 'AGATCGGAAGA'
##was used in this study
##qualityCheck: whether or not to filter for high-quality base reads ('yes' or 'no');
##'no' was used in this study
#output:
#1. tracking file reporting progress of script
#2. *_trimmed.fastq FASTQ file of filtered reads according to arguments above
#3. *_long.fastq FASTQ file of all reads that were too long to be considered further
#4. *_all.tab tab-separated file of all read sequences with adapters
#5. summary text file listing summary statistics related to filtering

import sys

import os

#read in FASTQ file

filePath = '/lustre/scratch/users/michael.nodine/seq/' #have to replace this with filePath of your system

folder = sys.argv[1]
lowerLimit = int(sys.argv[2].strip())
upperLimit = int(sys.argv[3].strip())
adapterSeq = sys.argv[4].strip()
qualityCheck = sys.argv[5].strip()

print adapterSeq

workDir = filePath + folder

outFile1 = workDir + '/tracking'
ofh1 = open(outFile1,'w')

fileList = os.listdir(workDir)

fileCount = 0

for file in fileList:
    if file.endswith('.fastq'):
        fastq = file
        fileCount += 1

print >>ofh1,'Currently processing fastq file through sRNA_step1.py......'

if fileCount > 1:
    print >>ofh1,'Error: more than one .fastq file detected'
else:
    print >>ofh1,'OK: only one .fastq file detected'

fastqFile = workDir + '/' + fastq

outFile = fastqFile[:-6] + '_trimmed.fastq'
longFile = fastqFile[:-6] + '_long.fastq'
allFile = fastqFile[:-6] + '_all.tab'

print >>ofh1,'fastqFile:',fastqFile
print >>ofh1,'outFile:',outFile
print >>ofh1,'allFile:',allFile

ofh = open(outFile,'w')
lfh = open(longFile,'w')
afh = open(allFile,'w')

totalTagNum = 0
a = 0
b = 4
noLinkerCounts = 0
LinkerCounts = 0
longEnoughCounts = 0
tooShortCounts = 0
goodQualCounts = 0
poorQualCounts = 0
tooLongCounts = 0
finalCount = 0
i = 0
markerCount = 0

#quality scroes based on Illumina 1.8+; poorQual_list is a list of quality scores below 30 (i.e an estimated probability of wrong base call = 0.001 [1/1000]) 

####################0###1###2###3###4###5###6###7###8###9###10##11##12##13##14##15##16##17##18##19##20##21##22##23##24##25##26##27##28##29##30##31##32##33##34##35##36##37##38##39##40

poorQual_list = ['!','"','#','$','%','&',"'",'(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>']#,'?','@','A','B','C','D']#,'E','F','G','H','I']

fastqFile = open(fastqFile)
lastLine = False
while lastLine == False:
    package = []
    totalTagNum = totalTagNum + 1
    for i in range(a, b):
        infileLine = fastqFile.readline()
        if infileLine == '':
            lastLine = True
        else:
            lastLine = False
            package.append(infileLine.rstrip())
        a = a + 4
        b = b + 4
    if lastLine != True:
        seq = package[1]
        qualScore = package[3]
        if adapterSeq != 'none':
            linkerMatch = seq.upper().find(adapterSeq)
            if linkerMatch == -1:
                noLinkerCounts += 1
            else:
                LinkerCounts += 1
                seq = seq[:linkerMatch].upper()
                print >>afh,seq,'\t',len(seq)
                #create a qualList and (optionally) search for characters within list for quality filter (note: in general this option is set to 'no' because filtering based on quality score is expected to increase # of false-negatives with little decrease of false-positives)
                qualScore = qualScore[:len(seq)]
                if len(seq) >= lowerLimit and 'N' not in seq:       #if 'N' is in 'trimmed' seq, then it cannot match to genome perfectly
                    longEnoughCounts += 1
                    if qualityCheck == 'yes':
                        qualList = list(qualScore)
                        qualCheck = 'Pass'
                        
                        for pQ in poorQual_list:
                            if pQ in qualList:
                                qualCheck = 'Fail'
                        if qualCheck == 'Pass':
                            goodQualCounts += 1
                            if len(seq) > upperLimit:
                                tooLongCounts += 1
                                print >>lfh, package[0]
                                print >>lfh, seq
                                print >>lfh, package[2]
                                print >>lfh, qualScore
                            else:
                                print >>ofh, package[0]
                                print >>ofh, seq
                                print >>ofh, package[2]
                                print >>ofh, qualScore
                                finalCount += 1
                        else:
                            poorQualCounts += 1
                    else:
                        if len(seq) > upperLimit:
                            tooLongCounts += 1
                            print >>lfh, package[0]
                            print >>lfh, seq
                            print >>lfh, package[2]
                            print >>lfh, qualScore
                        else:
                            print >>ofh, package[0]
                            print >>ofh, seq
                            print >>ofh, package[2]
                            print >>ofh, qualScore
                            finalCount += 1
                else:
                    tooShortCounts += 1

        else:
            seq = seq
            qualScore = qualScore
            
            if len(seq) >= lowerLimit and 'N' not in seq:
                longEnoughCounts += 1
                if qualityCheck == 'yes':
                    qualList = list(qualScore)
                    qualCheck = 'Pass'
                    
                    for pQ in poorQual_list:
                        if pQ in qualList:
                            qualCheck = 'Fail'
                    if qualCheck == 'Pass':
                        goodQualCounts += 1
                        if len(seq) > upperLimit:
                            tooLongCounts += 1
                            print >>lfh, package[0]
                            print >>lfh, seq
                            print >>lfh, package[2]
                            print >>lfh, qualScore
                        else:
                            print >>ofh, package[0]
                            print >>ofh, seq
                            print >>ofh, package[2]
                            print >>ofh, qualScore
                            finalCount += 1
                    else:
                        poorQualCounts += 1
                else:
                    if len(seq) > upperLimit:
                        tooLongCounts += 1
                        print >>lfh, package[0]
                        print >>lfh, seq
                        print >>lfh, package[2]
                        print >>lfh, qualScore
                    else:
                        print >>ofh, package[0]
                        print >>ofh, seq
                        print >>ofh, package[2]
                        print >>ofh, qualScore
                        finalCount += 1
            else:
                tooShortCounts += 1
            
fastqFile.close()
ofh.close()
lfh.close()
afh.close()

################################
outFile = workDir + '/summary'
################################
ofh = open(outFile,'w')

print >>ofh,'Total number of reads:','\t',totalTagNum - 1
print >>ofh,'Total number of reads with',adapterSeq,'of three prime adapter:','\t',LinkerCounts,'\t',round(LinkerCounts/float(totalTagNum - 1),2)
if adapterSeq != 'none':
    print >>ofh,'Total number of reads that are >=',lowerLimit,'\t',longEnoughCounts,'\t',round(longEnoughCounts/float(LinkerCounts),2),'\t',round(longEnoughCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads with quality scores >= 30:',goodQualCounts,'\t',round(goodQualCounts/float(longEnoughCounts),2),'\t',round(goodQualCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads with quality scores < 30:',poorQualCounts,'\t',round(poorQualCounts/float(longEnoughCounts),2),'\t',round(poorQualCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads that are >',upperLimit,'\t',tooLongCounts,'\t',round(tooLongCounts/float(LinkerCounts),2),'\t',round(tooLongCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of filtered reads:','\t',finalCount,'\t',round(finalCount/float(LinkerCounts),2),'\t',round(finalCount/float(totalTagNum - 1),2)
else:
    print >>ofh,'Total number of reads that are >=',lowerLimit,'\t',longEnoughCounts,'\t',round(longEnoughCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads with quality scores >= 30:',goodQualCounts,'\t',round(goodQualCounts/float(longEnoughCounts),2),'\t',round(goodQualCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads with quality scores < 30:',poorQualCounts,'\t',round(poorQualCounts/float(longEnoughCounts),2),'\t',round(poorQualCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of reads that are >',upperLimit,'\t',tooLongCounts,'\t',round(tooLongCounts/float(totalTagNum - 1),2)
    print >>ofh,'Total number of filtered reads:','\t',finalCount,'\t',round(finalCount/float(totalTagNum - 1),2)

ofh.close()

print >>ofh1,'sRNA_step1_var_qual.py complete'
ofh1.close()
