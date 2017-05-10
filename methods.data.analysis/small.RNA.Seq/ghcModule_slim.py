#ghcModule_slim.py
#Contains functions required for make_ghcsAndlDicts.py

class GenomicHit:
    # this may save some space by reducing the number of chromosome strings
    # that are associated with Locus instances (see __init__)
    __chrDict = dict()
    __strandDict = {'+':'+', '-':'-', '.':'.'}
    def __init__(self,name,seq,chromo,start,end,strand,readNum_raw,readNum_norm):
        self._name = name
        self._seq = seq
        coords = [start,end]
        #sort coords so that start < end
        coords.sort()
        if not(self.__chrDict.has_key(chromo)):
            self.__chrDict[chromo] = (chromo)
        self._chromo = self.__chrDict[chromo]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._strand = strand
        self._readNum_raw = int(readNum_raw)
        self._readNum_norm = readNum_norm
    def name(self):
        return self._name
    def seq(self):
        return self._seq
    def chromo(self):
        return self._chromo
    def start(self):
        return self._start
    def end(self):
        return self._end
    def strand(self):
        return self._strand
    def readNum_raw(self):
        return self._readNum_raw
    def readNum_norm(self):
        return self._readNum_norm
    def getAntisenseGenomicHit(self):
        if self._strand=='.':
            return self
        else:
            switch = {'+':'-', '-':'+'}
            name = self._name + '_AS'
            return GenomicHit(name,makeRC(self._seq),self._chromo,self._start,self._end,switch[self._strand],self._readNum_raw,self._readNum_norm)
    def overlaps(self,otherLocus):
        #are they on the same chromo? If not, then return False
        if self.chromo() != otherLocus.chromo():
            return False
        #do they match the ambidextrous or same strand? If not, return False
        elif not(self._strand == '.' or \
             otherLocus.strand() == '.' or \
             self.strand() == otherLocus.strand()): return False
        #do their start or end coordinates overlap? If not, return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end():
            return False
        else:
            return True
    #same as overlaps but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseGenomicHit().overlaps(otherLocus)


class GenomicHit_rpm:
    # this may save some space by reducing the number of chromosome strings
    # that are associated with Locus instances (see __init__)
    __chrDict = dict()
    __strandDict = {'+':'+', '-':'-', '.':'.'}
    def __init__(self,name,seq,chromo,start,end,strand,readNum_raw,readNum_normHit,readNum_normRpm):
        self._name = name
        self._seq = seq
        coords = [start,end]
        #sort coords so that start < end
        coords.sort()
        if not(self.__chrDict.has_key(chromo)):
            self.__chrDict[chromo] = (chromo)
        self._chromo = self.__chrDict[chromo]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._strand = strand
        self._readNum_raw = int(readNum_raw)
        self._readNum_normHit = readNum_normHit
        self._readNum_normRpm = readNum_normRpm
    def name(self):
        return self._name
    def seq(self):
        return self._seq
    def chromo(self):
        return self._chromo
    def start(self):
        return self._start
    def end(self):
        return self._end
    def strand(self):
        return self._strand
    def readNum_raw(self):
        return self._readNum_raw
    def readNum_normHit(self):
        return self._readNum_normHit
    def readNum_normRpm(self):
        return self._readNum_normRpm
    def getAntisenseGenomicHit(self):
        if self._strand=='.':
            return self
        else:
            switch = {'+':'-', '-':'+'}
            name = self._name + '_AS'
            return GenomicHit_rpm(name,makeRC(self._seq),self._chromo,self._start,self._end,switch[self._strand],self._readNum_raw,self._readNum_normHit,self._readNum_normRpm)
    def overlaps(self,otherLocus):
        #are they on the same chromo? If not, then return False
        if self.chromo() != otherLocus.chromo():
            return False
        #do they match the ambidextrous or same strand? If not, return False
        elif not(self._strand == '.' or \
             otherLocus.strand() == '.' or \
             self.strand() == otherLocus.strand()): return False
        #do their start or end coordinates overlap? If not, return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end():
            return False
        else:
            return True
    #same as overlaps but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseGenomicHit().overlaps(otherLocus)


class GenomicHitCollection:
    def __init__(self,gHits,windowSize):
        self.__chrToCoordToLoci = dict()
        self.__gHits = dict()
        self.__winSize = windowSize
        for gHit in gHits:
            self.__addHits(gHit)
    def __addHits(self,gHit):
        #check whether gHit is in ghc; if not then creates a chrKey (e.g. Ath_chr+) and adds to keyList
        #creates chrToCoordToLoci[chrKey] = ([range1]=(locusList1),[range2]=(locusList2))
        if not(self.__gHits.has_key(gHit)):
            self.__gHits[gHit] = None
            if gHit.strand()=='.':
                chrKeyList = [gHit.chromo()+'+', gHit.chromo()+'-']
            else:
                chrKeyList = [gHit.chromo()+gHit.strand()]
            #for each chrKey, checks whether key is already present in chrToCoordToLoci
            #if not, then it adds it
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)):
                    self.__chrToCoordToLoci[chrKey] = dict()
                #retrieves the range of coords
                #checks whether rangeKey exists in chrToCoordToLoci[chrKey] dict, if not then it adds it and
                #appends locus to corresponding list
                for n in self.__getKeyRange(gHit):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)):
                        self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(gHit)
    #helps with the getOverlap function
    def __subsetHelper(self,locus,strand):
        strand = strand.lower()
        matches = dict()
        senses = ['+','-']
        #depending on strand being searched, creates lambdas
        if locus.strand()=='.' or strand=='both':
            lamb = lambda s: True
        elif strand=='sense':
            lamb = lambda s: s==locus.strand()
        elif strand=='antisense':
            lamb = lambda s: s!=locus.strand()
        else: raise ValueError("sense value was inappropriate: '"+sense+"'.")
        #now search for appropriate chrToCoordToLoci dict entries; first depending on strand of locus
        #and search criteria
        for s in filter(lamb,senses):
            chrKey = locus.chromo()+s
            #if chrKey exists, then retrieve matches in matches dict
            if self.__chrToCoordToLoci.has_key(chrKey):
                #gets locus range and iterates through to find all potential hits within range
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        #if corresponding loci are within interval, then creates matches dict entries == None
                        for gHit in self.__chrToCoordToLoci[chrKey][n]:
                            matches[gHit] = None
        return matches.keys()
    #strand can be 'sense' (default), 'antisense' or 'both'
    #getOverlap returns all members of the collection that overlap the locus
    def getOverlap(self,locus,strand='sense'):
    #first use subsetHelper to identify all possible hits
        matches = self.__subsetHelper(locus,strand)
    #now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            #checks if potential matches identified with subset helper are real matches
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()#realMatches.keys()
    def __getKeyRange(self,gHit):
        start = gHit.start()/self.__winSize
        end = gHit.end()/self.__winSize + 1 ## add 1 because of the range
        return range(start,end)
    def get_readNum_raw(self):
##        ghc = self._ghc
        readNum_raw = 0
        for i in self.__gHits.keys():
            readNum_raw = readNum_raw + i.readNum_raw()
        return readNum_raw
    def get_readNum_norm(self):
##        ghc = self._ghc
        readNum_norm = 0
        for i in self.__gHits.keys():
            readNum_norm = readNum_norm + i.readNum_norm()
        return readNum_norm
    def get_hitNum(self):
##        ghc_list = self.get_ghc()
        hitNum = len(self.__gHits.keys())
        return hitNum
    def getHits(self):
        return self.__gHits.keys()
    def remove(self,old):
        if not(self.__gHits.has_key(old)): raise ValueError("requested locus isn't in collection")
        del self.__gHits[old]
        if old.strand()=='.': senseList = ['+','-']
        else: senseList = [old.strand()]
        for k in self.__getKeyRange(old):
            for sense in senseList:
                self.__chrToCoordToLoci[old.chromo()+sense][k].remove(old)

class GenomicHitCollection_rpm:
    def __init__(self,gHits,windowSize):
        self.__chrToCoordToLoci = dict()
        self.__gHits = dict()
        self.__winSize = windowSize
        for gHit in gHits:
            self.__addHits(gHit)
    def __addHits(self,gHit):
        #check whether gHit is in ghc; if not then creates a chrKey (e.g. Ath_chr+) and adds to keyList
        #creates chrToCoordToLoci[chrKey] = ([range1]=(locusList1),[range2]=(locusList2))
        if not(self.__gHits.has_key(gHit)):
            self.__gHits[gHit] = None
            if gHit.strand()=='.':
                chrKeyList = [gHit.chromo()+'+', gHit.chromo()+'-']
            else:
                chrKeyList = [gHit.chromo()+gHit.strand()]
            #for each chrKey, checks whether key is already present in chrToCoordToLoci
            #if not, then it adds it
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)):
                    self.__chrToCoordToLoci[chrKey] = dict()
                #retrieves the range of coords
                #checks whether rangeKey exists in chrToCoordToLoci[chrKey] dict, if not then it adds it and
                #appends locus to corresponding list
                for n in self.__getKeyRange(gHit):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)):
                        self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(gHit)
    #helps with the getOverlap function
    def __subsetHelper(self,locus,strand):
        strand = strand.lower()
        matches = dict()
        senses = ['+','-']
        #depending on strand being searched, creates lambdas
        if locus.strand()=='.' or strand=='both':
            lamb = lambda s: True
        elif strand=='sense':
            lamb = lambda s: s==locus.strand()
        elif strand=='antisense':
            lamb = lambda s: s!=locus.strand()
        else: raise ValueError("sense value was inappropriate: '"+sense+"'.")
        #now search for appropriate chrToCoordToLoci dict entries; first depending on strand of locus
        #and search criteria
        for s in filter(lamb,senses):
            chrKey = locus.chromo()+s
            #if chrKey exists, then retrieve matches in matches dict
            if self.__chrToCoordToLoci.has_key(chrKey):
                #gets locus range and iterates through to find all potential hits within range
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        #if corresponding loci are within interval, then creates matches dict entries == None
                        for gHit in self.__chrToCoordToLoci[chrKey][n]:
                            matches[gHit] = None
        return matches.keys()
    #strand can be 'sense' (default), 'antisense' or 'both'
    #getOverlap returns all members of the collection that overlap the locus
    def getOverlap(self,locus,strand='sense'):
    #first use subsetHelper to identify all possible hits
        matches = self.__subsetHelper(locus,strand)
    #now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            #checks if potential matches identified with subset helper are real matches
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()#realMatches.keys()

    #getContained returns all members of the collection that are contained by the locus
    def getContained(self,locus,strand='sense'):
        matches = self.__subsetHelper(locus,strand)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            for i in filter(lambda lcs: locus.contains(lcs), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: locus.containsAntisense(lcs), matches):
                realMatches[i] = None
        return realMatches.keys()

    def __getKeyRange(self,gHit):
        start = gHit.start()/self.__winSize
        end = gHit.end()/self.__winSize + 1 ## add 1 because of the range
        return range(start,end)
    def get_readNum_raw(self):
        readNum_raw = 0
        for i in self.__gHits.keys():
            readNum_raw = readNum_raw + i.readNum_raw()
        return readNum_raw
    def get_readNum_normHit(self):
        readNum_normHit = 0
        for i in self.__gHits.keys():
            readNum_normHit = readNum_normHit + i.readNum_normHit()
        return readNum_normHit
    def get_readNum_normRpm(self):
        readNum_normRpm = 0
        for i in self.__gHits.keys():
            readNum_normRpm = readNum_normRpm + i.readNum_normRpm()
        return readNum_normRpm
    def get_hitNum(self):
        hitNum = len(self.__gHits.keys())
        return hitNum
    def getHits(self):
        return self.__gHits.keys()
    def remove(self,old):
        if not(self.__gHits.has_key(old)): raise ValueError("requested locus isn't in collection")
        del self.__gHits[old]
        if old.strand()=='.': senseList = ['+','-']
        else: senseList = [old.strand()]
        for k in self.__getKeyRange(old):
            for sense in senseList:
                self.__chrToCoordToLoci[old.chromo()+sense][k].remove(old)
    def add(self,gHits,windowSize):
        #self.__chrToCoordToLoci = dict()
        #self.__gHits = dict()
        #self.__winSize = windowSize
        for gHit in gHits:
            self.__addHits(gHit)

class Locus:
    # this may save some space by reducing the number of chromosome strings
    # that are associated with Locus instances (see __init__)
    __chrDict = dict()
    __strandDict = {'+':'+', '-':'-', '.':'.'}
    def __init__(self,geneID,chromo,start,end,strand,gene_type,short_desc,long_desc,comp_desc):
        self._geneID = geneID
        coords = [start,end]
        #sort coords so that start < end
        coords.sort()
        if not(self.__chrDict.has_key(chromo)):
            self.__chrDict[chromo] = (chromo)
        self._chromo = self.__chrDict[chromo]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._strand = strand#self.__strandDict[strand]
        self._gene_type = gene_type
        self._short_desc = short_desc
        self._long_desc = long_desc
        self._comp_desc = comp_desc
    def geneID(self):
        return self._geneID
    def chromo(self):
        return self._chromo
    def start(self):
        return self._start
    def end(self):
        return self._end
    def strand(self):
        return self._strand
    def gene_type(self):
        return self._gene_type
    def short_desc(self):
        return self._short_desc
    def long_desc(self):
        return self._long_desc
    def comp_desc(self):
        return self._comp_desc
    def getAntisenseLocus(self):
        if self._strand=='.':
            return self
        else:
            switch = {'+':'-', '-':'+'}
            return Locus(self._geneID+'_AS',self._chromo,self._start,self._end,switch[self._strand],self._gene_type,self._short_desc,self._long_desc,self._comp_desc)
    #checks whether two loci overtlap
    def overlaps(self,otherLocus):
        #are they on the same chromo?
        if self.chromo() != otherLocus.chromo():
            return False
        #do they match ambidextrous or the same strand? If not, return False
        elif not(self._strand == '.' or \
             otherLocus.strand() == '.' or \
             self.strand() == otherLocus.strand()): return False
        #do they overlap
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end():
            return False
        else:
            return True
    #checks whether all the nts of given locus overlap with self locus
    def contains(self,otherLocus):
        if self.chromo()!=otherLocus.chromo(): return False
        elif not(self._strand=='.' or \
                 otherLocus.strand()=='.' or \
                 self.strand()==otherLocus.strand()): return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end(): return False
        else: return True
    # same as contains, but considers the opposite strand
    def containsAntisense(self,otherLocus):
        return self.getAntisenseLocus().contains(otherLocus)
    #same as overlaps but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseLocus().overlaps(otherLocus)


class Locus_miRNA:
    __chrDict = dict()
    __strandDict = {'+':'+', '-':'-', '.':'.'}
    def __init__(self,geneID,chromo,start,end,strand,gene_type,common,AGI_ID,MB_ID,family,short_desc,matList,pTargets,vTargets):
        self._geneID = geneID
        coords = [start,end]
        #sort coords so that start < end
        coords.sort()
        if not(self.__chrDict.has_key(chromo)):
            self.__chrDict[chromo] = (chromo)
        self._chromo = self.__chrDict[chromo]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._strand = strand#self.__strandDict[strand]
        self._gene_type = gene_type
        self._common = common
        self._AGI_ID = AGI_ID
        self._MB_ID = MB_ID
        self._family = family
        self._short_desc = short_desc
        self._matList = matList
        self._pTargets = pTargets
        self._vTargets = vTargets
    def geneID(self):
        return self._geneID
    def chromo(self):
        return self._chromo
    def start(self):
        return self._start
    def end(self):
        return self._end
    def strand(self):
        return self._strand
    def gene_type(self):
        return self._gene_type
    def common(self):
        return self._common
    def AGI_ID(self):
        return self._AGI_ID
    def MB_ID(self):
        return self._MB_ID
    def family(self):
        return self._family
    def short_desc(self):
        return self._short_desc
    def matList(self):
        return self._matList
    def pTargets(self):
        return self._pTargets
    def vTargets(self):
        return self._vTargets
    def getAntisenseLocus(self):
        if self._strand=='.':
            return self
        else:
            switch = {'+':'-', '-':'+'}
            return Locus_miRNA(self._geneID+'_AS',self._chromo,self._start,self._end,switch[self._strand],self._gene_type,self._common,self._AGI_ID,self._MB_ID,self._family,self._short_desc,self._matList,self._pTargets,self._vTargets)
    #checks whether two loci overtlap
    def overlaps(self,otherLocus):
        #are they on the same chromo?
        if self.chromo() != otherLocus.chromo():
            return False
        #do they match ambidextrous or the same strand? If not, return False
        elif not(self._strand == '.' or \
             otherLocus.strand() == '.' or \
             self.strand() == otherLocus.strand()): return False
        #do they overlap
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end():
            return False
        else:
            return True
    #checks whether all the nts of given locus overlap with self locus
    def contains(self,otherLocus):
        if self.chromo()!=otherLocus.chromo(): return False
        elif not(self._strand=='.' or \
                 otherLocus.strand()=='.' or \
                 self.strand()==otherLocus.strand()): return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end(): return False
        else: return True
    # same as contains, but considers the opposite strand
    def containsAntisense(self,otherLocus):
        return self.getAntisenseLocus().contains(otherLocus)
    #same as overlaps but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseLocus().overlaps(otherLocus)


class LocusCollection:
    def __init__(self,loci,windowSize):
        self.__chrToCoordToLoci = dict()
        self.__loci = dict()
        self.__winSize = windowSize
        for lcs in loci:
            self.__addLocus(lcs)
    def __addLocus(self,lcs):
        #check whether loci is in locus collection; if not then creates a chrKey (e.g. Ath_chr+) and adds to keyList
        #creates chrToCoordToLoci[chrKey] = ([range1]=(locusList1),[range2]=(locusList2))
        if not(self.__loci.has_key(lcs)):
            self.__loci[lcs] = None
            if lcs.strand()=='.':
                chrKeyList = [lcs.chromo()+'+', lcs.chromo()+'-']
            else:
                chrKeyList = [lcs.chromo()+lcs.strand()]
            #for each chrKey, checks whether key is already present in chrToCoordToLoci
            #if not, then it adds it
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)):
                    self.__chrToCoordToLoci[chrKey] = dict()
                #retrieves the range of coords
                #checks whether rangeKey exists in chrToCoordToLoci[chrKey] dict, if not then it adds it and
                #appends locus to corresponding list
                for n in self.__getKeyRange(lcs):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)):
                        self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(lcs)

    #helps with the getOverlap function
    def __subsetHelper(self,locus,strand):
        matches = dict()
        senses = ['+','-']
        #depending on strand being searched, creates lambdas
        if locus.strand()=='.' or strand=='both':
            lamb = lambda s: True
        elif strand=='sense':
            lamb = lambda s: s==locus.strand()
        elif strand=='antisense':
            lamb = lambda s: s!=locus.strand()
        #now search for appropriate chrToCoordToLoci dict entries; first depending on strand of locus
        #and search criteria
        for s in filter(lamb,senses):
            chrKey = locus.chromo()+s
            #if chrKey exists, then retrieve matches in matches dict
            if self.__chrToCoordToLoci.has_key(chrKey):
                #gets locus range and iterates through to find all potential hits within range
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        #if corresponding loci are within interval, then creates matches dict entries == None
                        for lcs in self.__chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None
        return matches.keys()

    #strand can be 'sense' (default), 'antisense' or 'both'
    #getOverlap returns all members of the collection that overlap the locus
    def getOverlap(self,locus,strand='sense'):
    #first use subsetHelper to identify all possible hits
        matches = self.__subsetHelper(locus,strand)
    #now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            #checks if potential matches identified with subset helper are real matches
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()
    #getContained returns all members of the collection that are contained by the locus
    def getContained(self,locus,strand='sense'):
        matches = self.__subsetHelper(locus,strand)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            for i in filter(lambda lcs: locus.contains(lcs), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: locus.containsAntisense(lcs), matches):
                realMatches[i] = None
        return realMatches.keys()
    # returns all members of the collection that contain the locus
    def getContainers(self,locus,strand='sense'):
        matches = self.__subsetHelper(locus,strand)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if strand=='sense' or strand=='both':
            for i in filter(lambda lcs: lcs.contains(locus), matches):
                realMatches[i] = None
        if strand=='antisense' or strand=='both':
            for i in filter(lambda lcs: lcs.containsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()
    def __getKeyRange(self,locus):
        start = locus.start()/self.__winSize
        end = locus.end()/self.__winSize + 1 ## add 1 because of the range
        return range(start,end)
    def getLoci(self):
        return self.__loci.keys()
    def getGeneIDs(self):
        idList = []
        lc_list = self.getLoci()
        for lcs in lc_list:
            if lcs.geneID() not in idList:
                idList.append(lcs.geneID())
        return idList

#makeRC converts a seq into its reverse complement

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
