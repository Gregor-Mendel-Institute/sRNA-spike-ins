#methods.sRNA.spike.in.design.step2.endogenous.miRNAs.py
#
#input: mature.miRNA.seqs.top50percent.fa
#output: 1) folding structures calculated by RNAfold, 2) minimum free energies (MFEs) calculated by RNAfold

import sys, os

inFile = 'mature.miRNA.seqs.top50percent.fa'
##ifh = open(inFile)
##inLines = ifh.readlines()
##ifh.close()
##
##faDict = {}
##
##for i in range(0,len(inLines),2):
##    name = inLines[i].strip()
##    name = name.lstrip('>')
##    seq = inLines[i+1].strip()
##    faDict[name] = (seq)


cmd = 'cat ' + inFile + '| RNAfold -T 4 --noPS > ' + 'mature.miRNA.seqs.top50percent_folded'

os.system(cmd)


#retrieve mfes

outFile = 'mature.miRNA.seqs.top50percent_mfes'
ofh = open(outFile,'w')

inFile = 'mature.miRNA.seqs.top50percent_folded'
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
    print >>ofh,mfe

ofh.close()    

