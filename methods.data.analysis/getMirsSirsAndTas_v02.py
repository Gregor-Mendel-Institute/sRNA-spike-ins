#getMirsSirsAndTas_v02.py
#extracts levels (in RPM) for miRNA families, tasiRNA families, and 20-22 nt or 23-24 nt reads contained within annotated TEs (siRNAs)the following and prints these to seperate files for further processing:

outDir = "data.for.graphs/"

sampleList = ['col0_leaf1','col0_leaf2','col0_fb1','col0_fb2','d234_fb1','d234_fb2']

teList = []

teFile1 = 'Ath_annotations/te_AGIs'
teFile2 = 'Ath_annotations/teg_AGIs'

ifh1 = open(teFile1)
ifh2 = open(teFile2)

inLines1 = ifh1.readlines()
inLines2 = ifh2.readlines()

ifh1.close()
ifh2.close()

teList_anno = []

for line in inLines1[2:]:
    cols = line.split('\t')
    te = cols[0].strip()
    if te not in teList:
        teList_anno.append(te)

for line in inLines2[2:]:
    cols = line.split('\t')
    te = cols[0].strip()
    if te not in teList:
        teList_anno.append(te)

print 'Length of teList_anno:',len(teList_anno)


for sample in sampleList:
    print sample
    miRfile = sample + '/miRNA/families_byRpm'
    tasiRfile = sample + '/tasiRNA/families_byRpm'
    siRfile = sample + '/gHits_singleGeneTypes/te_gHits'

    outFile_miR = outDir + sample + '_miR_fams'
    outFile_tasiR = outDir + sample + '_tasiR_fams'
    outFile_siR_1 = outDir + sample + '_siR_24'
    outFile_siR_2 = outDir + sample + '_siR_21'

    ofh_miR = open(outFile_miR,'w')
    ofh_tasiR = open(outFile_tasiR,'w')
    ofh_siR_1 = open(outFile_siR_1,'w')
    ofh_siR_2 = open(outFile_siR_2,'w')

    print >>ofh_miR,'family','\t','rpm'
    print >>ofh_tasiR,'family','\t','rpm'
    print >>ofh_siR_1,'TE','\t','rpm'
    print >>ofh_siR_2,'TE','\t','rpm'

    ifh_miR = open(miRfile)
    inLines_miR = ifh_miR.readlines()
    ifh_miR.close()

    for line in inLines_miR[1:]:
        cols = line.strip().split('\t')
        family = cols[0].strip()
        rpm = cols[1].strip()
        print >>ofh_miR,family,'\t',rpm
            
    ifh_tasiR = open(tasiRfile)
    inLines_tasiR = ifh_tasiR.readlines()
    ifh_tasiR.close()

    for line in inLines_tasiR[1:]:
        cols = line.strip().split('\t')
        family = cols[0].strip()
        rpm = cols[1].strip()
        print >>ofh_tasiR,family,'\t',rpm

    ifh_te = open(siRfile)
    inLines_te = ifh_te.readlines()

    te_21 = {}      #[TE] = (rpm)
    te_24 = {}

    #prepopulate te dicts with known TEs
    for te in teList_anno:
        te_21[te] = (0)
        te_24[te] = (0)

    for line in inLines_te[1:]:
        cols = line.strip().split('\t')
        seq = cols[1].strip()
        if len(seq) in range(20,23):
            rpm = float(cols[8].strip())
            teList = cols[9].strip().split('$')
            teList_new = []
            for te in teList[1:]:
                teList_new.append(te.split('_')[0])

            for te in teList_new:
                rpm_new = te_21[te] + rpm
                te_21[te] = (rpm_new)
        if len(seq) in range(23,25):
            rpm = float(cols[8].strip())
            teList = cols[9].strip().split('$')
            teList_new = []
            for te in teList[1:]:
                teList_new.append(te.split('_')[0])

            for te in teList_new:
                rpm_new = te_24[te] + rpm
                te_24[te] = (rpm_new)

    print 'Number of TEs with 20-22 nt siRNAs:',len(te_21.keys())
    print 'Number of TEs with 23-24 nt siRNAs:',len(te_24.keys())
    
    for te in te_21.keys():
        print >>ofh_siR_1,te,'\t',te_21[te]
        
    for te in te_24.keys():
        print >>ofh_siR_2,te,'\t',te_24[te]
            
    ofh_miR.close()
    ofh_tasiR.close()
    ofh_siR_1.close()
    ofh_siR_2.close()
