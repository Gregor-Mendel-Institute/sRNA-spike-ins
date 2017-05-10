#getMirsSirsAndTas_v03.py
#extracts levels (in RPM) for the following and prints to corresponding outfiles:
#   miRNAS: 20-22 nt reads that are contained within TAIR10/miRBase21 annotated miRNAs
#   tasiRNAS: 20-22 nt reads that are contained within TAIR10/automated annotated tasiRNAs
#   20-22 nt siRNAs: 20-22 nt reads that are contained within TE-related fragments
#   23-24 nt siRNAs: 23-24 nt reads that are contained within TE-related fragments

outDir = "data.for.graphs/"

sampleList = ['col0_leaf1','col0_leaf2','col0_fb1','col0_fb2','d234_fb1','d234_fb2']

for sample in sampleList:
    print sample
    miRfile = sample + '/miRNA/matureFound'
    tasiRfile = sample + '/tasiRNA/matureFound'
    siRfile = sample + '/gHits_all/te_gHits'

    outFile_miR = outDir + sample + '_miR_mature'
    outFile_tasiR = outDir + sample + '_tasiR_mature'
    outFile_siR_1 = outDir + sample + '_siR_24_v2'
    outFile_siR_2 = outDir + sample + '_siR_21_v2'

    ofh_miR = open(outFile_miR,'w')
    ofh_tasiR = open(outFile_tasiR,'w')
    ofh_siR_1 = open(outFile_siR_1,'w')
    ofh_siR_2 = open(outFile_siR_2,'w')

    print >>ofh_miR,'rpm'
    print >>ofh_tasiR,'rpm'
    print >>ofh_siR_1,'rpm'
    print >>ofh_siR_2,'rpm'

    ifh_miR = open(miRfile)
    inLines_miR = ifh_miR.readlines()
    ifh_miR.close()

    for line in inLines_miR[1:]:
        cols = line.strip().split('\t')
        seq = cols[5].strip()
        rpm = cols[8].strip()
        if len(seq) in range(20,23):
            print >>ofh_miR,rpm
            
    ifh_tasiR = open(tasiRfile)
    inLines_tasiR = ifh_tasiR.readlines()
    ifh_tasiR.close()

    for line in inLines_tasiR[1:]:
        cols = line.strip().split('\t')
        seq = cols[5].strip()
        rpm = cols[8].strip()
        if len(seq) in range(20,23):
            print >>ofh_tasiR,rpm
            
    ifh_siR = open(siRfile)
    inLines_siR = ifh_siR.readlines()
    ifh_siR.close()

    for line in inLines_siR[1:]:
        cols = line.strip().split('\t')
        seq = cols[1].strip()
        rpm = cols[8].strip()
        if len(seq) in range(23,25):
            print >>ofh_siR_1,rpm
        elif len(seq) in range(20,23):
            print >>ofh_siR_2,rpm
    
    ofh_miR.close()
    ofh_tasiR.close()
    ofh_siR_1.close()
    ofh_siR_2.close()
