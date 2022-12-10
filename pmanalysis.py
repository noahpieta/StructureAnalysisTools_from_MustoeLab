

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
        
        self.sourcename = corrfile

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
    

    def ppvsens_duplex(self, ctobj, ptype=1, exact=False, profile=None, mask=False, verbal=False, printFP=False):
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
        
        mask = True will exclude masked regions in ctobj from ppv/sens calcs

        """
        
        if mask:
            print(len(ctobj.mask), len(ctobj.mask)-sum(ctobj.mask))

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
        

        predpairs = []
        
        for c in corrlist:
            
            pair = (c[0], c[1])
            if not (mask and masked(pair, ctobj.mask, self.window)):
                predpairs.append(pair)
            elif verbal:
                print('Skipped {0}'.format(pair))

        predpairs = set(predpairs)

        
        # for sensitivity calculation, get list of helices
        helices = ctobj.extractHelices(fillPairs=False)
        
        # need to shift helix indices to match window-shifting in PairMap
        shifted_helices = {}
        knowndup = []
        for h,pairs in helices.items():
            
            temphelix = []
            
            for p in pairs[:-(self.window-1)]:
                
                shiftedpair = (p[0], p[1]-self.window+1)
                
                # skip if masked
                if mask and masked(shiftedpair, ctobj.mask, self.window):
                    continue

                if profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)

            
            # this deals with helices smaller than the window
            # for calculations, add in "dummies" that starts and end 1 bp upstream/downstream
            if len(pairs) < self.window:
                
                shiftedpair = (pairs[0][0]-1, pairs[0][1]-self.window+2)
                # skip if masked
                if mask and masked(shiftedpair, ctobj.mask, self.window):
                    pass
                elif profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)

                shiftedpair = (pairs[0][0], pairs[0][1]-self.window+1)
                if mask and masked(shiftedpair, ctobj.mask, self.window):
                    pass
                elif profile is None or hasdata(shiftedpair, profile, self.window):
                    temphelix.append(shiftedpair)
            

            if len(temphelix) > 0:
                knowndup.extend( temphelix )
                shifted_helices[h] = temphelix
                if verbal:
                    print(h, temphelix)
            elif verbal:
                print(h, 'skipped', pairs)
            #    print "WARNING: Helix starting at pair {0} skipped because of no data".format(pairs[0])

        
        knowndup = set(knowndup)

        # compute sens
        senset = set()
        for h in shifted_helices:
            for c in shifted_helices[h]:
                if c in predpairs or nonexactMatch(c, predpairs, allowedoffset):
                    senset.add(h)
                    break

        # compute ppv
        ppvset = set()
        for c in sorted(predpairs):
            if c in knowndup or nonexactMatch(c, knowndup, allowedoffset):
                ppvset.add(c)
            elif printFP:
                print('FP: {}'.format(c))    
        

        if len(predpairs)==0:
            ppv = 0.0
        else:
            ppv = float(len(ppvset))/len(predpairs)
        
        if len(shifted_helices)==0:
            sens = 0.0
        else:
            sens = float(len(senset))/len(shifted_helices)

        return ppv, sens
        

                   


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


def masked(pair, mask, window):
    
    s1 = sum(mask[pair[0]-1:pair[0]-1+window])
    s2 = sum(mask[pair[1]-1:pair[1]-1+window])
        
    if max(s1,s2) > window/2:
        return True
    else:
        return False
        
        
    



if __name__ == '__main__':
    
    import argparse
    import RNAtools2 as RNAtools

    prs = argparse.ArgumentParser()
    
    prs.add_argument('pmfile', help='path of pair mapper file')
    prs.add_argument('ctfile', help='path of reference ct file')
    prs.add_argument('--ptype', default=1, type=int, help='Which correlations to compute Sens/PPV for. Options are 1 (primary only;default), 2 (secondary only), 3 (primary and secondary), 0 (everything; not recommended)')
    prs.add_argument('--dms', help='path of dms reactivity file')
    prs.add_argument('--mask', action='store_true',help='Mask out CT-specified regions when computing ppv/sens')
    prs.add_argument('--verbal', action='store_true',help='Print each correlation and its status')
        
    args = prs.parse_args()


    pm = PairMap(args.pmfile)   
    ct = RNAtools.CT(args.ctfile, filterNC=True, filterSingle=True)
    
    if args.dms:
        profile,seq = RNAtools.readSHAPE(args.dms)
        if len(profile)!=len(ct.ct):
            raise IndexError('Profile file not the same length as the CT!')
    else:
        profile = None
    
    
    p,s = pm.ppvsens_duplex(ct, ptype=args.ptype, exact=False, profile=profile, mask=args.mask, verbal=args.verbal)
    print("PPV={0:.0f}  Sens={1:.0f}".format(p*100, s*100))

    #print pm.ppvsens_duplex(ct, ptype=1, exact=True, profile=profile)




