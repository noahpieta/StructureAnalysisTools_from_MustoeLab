

#######################################################################
#
# PairMaP object and related tools for post-analysis of PAIR-MaP data
#
# Anthony Mustoe
# Copyright 2018
# 
# This file is licensed under the terms of the MIT license 
#
# Change-log
#
#
#######################################################################




class PairMap(object):
    """This class is a container for correlations yielded by a RING experiment"""

    def __init__(self, corrfile=None):
        
        if corrfile is not None:
            self.readPairFile(corrfile)

    def readPairFile(self, corrfile):
        
        self.primary = []
        self.secondary = []
        self.remainder = []

        with open(corrfile) as inp:
            
            header = inp.readline().split()
            self.molsize = int(header[0])
            self.window = int(header[1].split('=')[1])
      
            # pop off the second header line
            inp.readline()
            
            for line in inp:
                spl = line.split()

                ptype = int(spl[3])
                c = (int(spl[0]), int(spl[1]), float(spl[2]), float(spl[4]))

                if ptype == 1:
                    self.primary.append(c)
                elif ptype == 2:
                    self.secondary.append(c)
                else:
                    self.remainder.append(c)
        
    
  
    #def ppvsens(self, corrlist):
    #    
    #    exit('ppvsens not working')
    #
    #    if self.knownstructure is None:
    #        raise AttributeError('knownstructure has not been initialized!')
    #    
    #    # calculate how many non-repetitive pairs match known pairs
    #    totalknown = 0.0
    #    for c in corrs_nr:
    #        if self.knownstructure.ct[c[0]-1] == c[1]:
    #            totalknown += 1
    #
    #    knownpairs = len(self.knownstructure.pairList())
    #    
    #    sens = totalknown / knownpairs
    #    ppv = totalknown / len(corrs_nr)    
    #
    #    return ppv, sens
    #
    

    def ppvsens_duplex(self, ctobj, ptype=1, exact=False, profile=None):
        """Compute ppv&sens from a duplex perspective relative to reference ct

        ptype = 1/2/3/0
                1 -> only primary
                2 -> only secondary
                3 -> 1+2
                0 -> everything (including remainder)

        exact = True/False
                if False, allow duplex to be shifted +/-1 from correct register 

        profile = if provided, reactivity profile is used to define regions with 
                  no data that are then excluded from ppv/sens calcs
        """
        

        if exact:
            allowedoffset = 0
        else:
            allowedoffset = 1


        corrlist = []
        if ptype == 1:
            corrlist = self.primary
        elif ptype == 2:
            corrlist = self.secondary
        elif ptype == 3:
            corrlist = self.primary + self.secondary
        elif ptype == 0:
            corrlist = self.primary + self.secondary + self.remainder
        else:
            raise ValueError('ptype={} is not supported. Please select 1/2/3/0'.format(ptype))
        
        predpairs = set([(c[0],c[1]) for c in corrlist])

        
        # for sensitivity calculation, get list of helices
        helices = ctobj.extractHelices(fillPairs=False)
        
        # need to shift helix indices to match window-shifting in PairMap
        shifted_helices = {}
        knowndup = []
        

        for h,pairs in helices.items():
            
            temphelix = []
            
            for p in pairs[:-(self.window-1)]:
                
                shiftedpair = (p[0], p[1]-self.window+1)
                
                if profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)

            
            # this deals with helices smaller than the window
            # for calculations, add in "dummies" that starts and end 1 bp upstream/downstream
            if len(pairs) < self.window:
                
                #print 'Warning! {0}-bp helix included in calculation; pairs={1}'.format(len(helices[h]), helices[h])

                shiftedpair = (pairs[0][0]-1, pairs[0][1]-self.window+2)
                if profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)

                shiftedpair = (pairs[0][0], pairs[0][1]-self.window+1)
                if profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)

            
            if len(temphelix) > 0:
                knowndup.extend( temphelix )
                shifted_helices[h] = temphelix
            #else:
            #    print "helix starting at pair {0} skipped because of no data".format(pairs[0])


        knowndup = set(knowndup)

        # compute sens
        sens = set()
        for h in shifted_helices:
            for c in shifted_helices[h]:
                if c in predpairs or nonexactMatch(c, predpairs, allowedoffset):
                    sens.add(h)
                    break

        # compute ppv
        ppv = set()
        for c in sorted(predpairs):
            if c in knowndup or nonexactMatch(c, knowndup, allowedoffset):
                ppv.add(c)
        
        #print len(ppv), len(predpairs), len(sens), len(shifted_helices)

        return float(len(ppv))/len(predpairs), float(len(sens))/len(shifted_helices)
        

                   


def nonexactMatch(corr, corrlist, allowedoffset=0):
    
    for x in range(-allowedoffset, allowedoffset+1):
        for y in range(-allowedoffset, allowedoffset+1):
            if (corr[0]+x,corr[1]+y) in corrlist:
                return True
    return False




def hasdata(pair, profile, window):

    missing = [0,0]
    for i in range(window):
        
        r0 = profile[pair[0]-1+i] # pair is 1-indexed, profile isn't
        r1 = profile[pair[1]-1+i] 
        
        if r0!=r0 or r0<-10:
            missing[0] += 1
        if r1!=r1 or r1<-10:
            missing[1] += 1
    
    if max(missing) == window:
        return False
    else:
        return True



    

if __name__ == '__main__':
    
    import argparse
    import RNAtools2 as RNAtools

    prs = argparse.ArgumentParser()
    
    prs.add_argument('pmfile', help='path of pair mapper file')
    prs.add_argument('ctfile', help='path of reference ct file')
    prs.add_argument('--dms', help='path of dms reactivity file')
        
    args = prs.parse_args()


    pm = PairMap(args.pmfile)   
    ct = RNAtools.CT(args.ctfile, filterNC=True)


    if args.dms:
        profile,seq = RNAtools.readSHAPE(args.dms)
    else:
        profile = None
    
    
    print pm.ppvsens_duplex(ct, ptype=1, exact=False, profile=profile)
    #print pm.ppvsens_duplex(ct, ptype=1, exact=True, profile=profile)




