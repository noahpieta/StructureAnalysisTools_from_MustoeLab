#!/usr/bin/python

# ---------------------------------------------------------------------------------------
# Flexible Arc Plotting code
# Anthony Mustoe
# Weeks lab, UNC
# 2016, 2017, 2018
# 
# Modifications
# Version 2.1 with plot profile option
#
# ---------------------------------------------------------------------------------------


import sys, os, math, argparse

import RNAtools2 as RNAtools

import matplotlib 
#matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['xtick.minor.size'] = 4 
matplotlib.rcParams['xtick.minor.width'] = 1
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Arial'


import matplotlib.pyplot as plot
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np


class ArcLegend(object):
    """Container for Legend for arc plot"""
    
    def __init__(self, title=None, colors=[], labels=[], msg=None):
        
        assert len(colors)==len(labels), "colors and labels must be the same size" 

        self.title=title
        self.colors = colors
        self.labels = labels
        self.msg = msg
    
    
    def append(self, t, colors, labels):
        
        if self.title is not None and t is not None:
            self.title+=' & '+t
        elif self.title is None:
            self.title = t

        self.colors.extend(colors)
        self.labels.extend(labels)


    def drawLegend(self, ax, xbound, ybound):
        """Draw legend on axes, computing location based on axes bounds
        ax     = axes object to draw on
        xbound = xboundary of axis
        ybound = yboundary of axis
        """
        
        spacing = 3
        
        # determine bounding box
        nlines = len(self.colors) + int(self.title is not None) + 1
        height = spacing*nlines       
        width = 15
        if self.title is not None and len(self.title) > width:
            width = len(self.title)
        
        #xloc = xbound-width
        #if xloc < 0:
        #    xloc = 0
        #if ybound < 0:
        #    yloc = ybound+height+2
        #    ax.add_patch(patches.Rectangle( (xloc-2, yloc-height+4), width, height, fc='white', lw=1.0))
        #else:
        #    yloc = ybound-height+2
        #    ax.add_patch(patches.Rectangle( (xloc-2, yloc+4), width, height, fc='white', lw=1.0))
        
        xloc = xbound
        yloc = ybound
        

        # write the title
        if self.title is not None:
            ax.text(xloc, yloc, self.title, horizontalalignment='left',
                    size="11",weight="medium", color="black")
            yloc -= spacing


        for i, c in enumerate(self.colors):

            # self.colors is a list of dict patch specifications, so c is a dict
            ax.add_patch(patches.Rectangle( (xloc, yloc), 1.5, 1.5, color=c, clip_on=False ) )

            # now add label
            ax.text(xloc+2.5, yloc, self.labels[i], horizontalalignment='left',
                    size="8",weight="normal", color="black")
            yloc -= spacing
    
        if self.msg:
            ax.text(xloc, yloc, self.msg, horizontalalignment='left',
                    size="8",weight="normal", color="black")



class ArcPlot(object):

    def __init__(self, title='', fasta=None):
        
        self.title = title

        self.seq = ''
        if fasta:
            self.addFasta(fasta)


        self.seqcolors = ['black']
        self.drawseq = True
        self.grid = False

        self.handleLength = 4*(math.sqrt(2)-1)/3
        self.adjust = 1.4

        self.topPatches = []
        self.botPatches = []
        self.height = [0, 0] # the figure height

        self.intdistance = None

        self.reactprofile = None
        self.reactprofileType = 'SHAPE'

        self.toplegend = None
        self.botlegend = None


    def addFasta(self, fastafile):
        
        read = False
        with open(fastafile) as inp:
            for line in inp:
   
                if line[0]=='>':
                    read = True
                    continue

                line = line.strip()
    
                if read and len(line) > 0:
                    self.seq += line

                                     



    def addArcPath(self, outerPair, innerPair=None, panel=1, color = 'black', alpha=0.5, window=1):
        """ add the arPath object for a given set of parameters
        If only outerPair is passed, innerPair is computed from window
        If both outerPair and innerPair are passed, window is not used
        """
        
        if innerPair is None:
            innerPair = [outerPair[0]+0.5 + window-1, outerPair[1]-0.5]
            outerPair = [outerPair[0]-0.5, outerPair[1]+0.5 + window-1]
        else:
            innerPair = [innerPair[0]+0.5, innerPair[1]-0.5]
            outerPair = [outerPair[0]-0.5, outerPair[1]+0.5]
            
        innerRadius = (innerPair[1] - innerPair[0])/2.0
        outerRadius = (outerPair[1] - outerPair[0])/2.0
        
        verts = []
        codes = []

        # outer left
        verts.append( (outerPair[0], 0) )
        codes.append( Path.MOVETO )

        # outer left control 1
        verts.append( (outerPair[0], panel*self.handleLength*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer left control 2
        verts.append( (outerPair[0]+outerRadius*(1-self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer center
        verts.append( (outerPair[0]+outerRadius, panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 1
        verts.append( (outerPair[0]+outerRadius*(1+self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 2
        verts.append( (outerPair[1], panel*self.handleLength*outerRadius) )
        codes.append( Path.LINETO )
                       
        # outer right
        verts.append( (outerPair[1], 0) )
        codes.append( Path.LINETO )
                        
        # inner right
        verts.append( (innerPair[1], 0) )
        codes.append( Path.LINETO )
            
        # inner right control 1
        verts.append( (innerPair[1], panel*self.handleLength*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner right control 2
        verts.append( (innerPair[0]+innerRadius*(1+self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner center
        verts.append( (innerPair[0]+innerRadius, panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 1 
        verts.append( (innerPair[0]+innerRadius*(1-self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 2
        verts.append( (innerPair[0], panel*self.handleLength*innerRadius) )
        codes.append( Path.LINETO )
 
        # inner left
        verts.append( (innerPair[0], 0) )
        codes.append( Path.LINETO )
 
        # outer left duplicate 
        verts.append( (outerPair[0], 0) )
        codes.append( Path.CLOSEPOLY )
 
        # rescale verts

        if panel == 1:
            adval = self.adjust
        else:
            adval = -(self.adjust-1)
        
        # move arcs away from x-axis
        verts = [ (x,y+adval) for x,y in verts ]
         

        indpath = Path(verts, codes)
        patch = patches.PathPatch(indpath, facecolor=color, alpha=alpha,
                                  linewidth=0, edgecolor='none')
        
        if panel == 1:
            self.topPatches.append(patch)
            if outerRadius > self.height[1]:
                self.height[1] = outerRadius
        
        else:
            self.botPatches.append(patch)
            if outerRadius > self.height[0]:
                self.height[0] = outerRadius



    def plotProfile(self, ax, bounds=None, colthresh = (-10, 0.4, 0.85, 3), heightscale=None):
        """Add a reactivity profile to the axes (self.reactprofile needs to be set)
        
        ax          = axes object to add plot. Expected to be top axis
        colthresh   = thresholds to bin reactivity data. 
                      First element indicates "No-data lower bound"
                      Last element indicates "Max reactivity" -- values are trucated above this
        heightscale = scaling to control height of bars
        """
        
        if self.reactprofile is None:
            return
        
        # check to see if DMS and colthresh defaults not overridden
        if self.reactprofileType == 'DMS' and colthresh == (-10,0.4,0.85,3):
            colthresh = (-10, 0.15, 0.4, 1)
            

        xvals = [ [] for i in range(4) ]
        yvals = [ [] for i in range(4) ]
        
        if heightscale is None:
            heightscale = max(4, min(10, len(self.reactprofile)/50.))
            heightscale = min(max(self.height)/4., heightscale)
            
            if self.reactprofileType == 'DMS': # adjust for compressed ploting range
                heightscale *= 2


        for x,y in enumerate(self.reactprofile):

            if bounds is not None and x<bounds[0] or x>bounds[1]:
                continue

            if y is None or y != y or y<colthresh[0]:
                xvals[0].append(x+1)
                yvals[0].append(0.6*heightscale+self.adjust)
            elif y < colthresh[1]:
                xvals[1].append(x+1)
                yvals[1].append(y*heightscale)
            elif y < colthresh[2]:
                xvals[2].append(x+1)
                yvals[2].append(y*heightscale)
            else:
                xvals[3].append(x+1)
                if y > colthresh[3]:
                    yvals[3].append(colthresh[3]*heightscale)
                else:
                    yvals[3].append(y*heightscale)

        
        ax.bar(xvals[0], yvals[0], alpha=0.7, linewidth=0, color=(179./255, 171./255, 148./255),
               align='center', clip_on=False, bottom=-0.3*heightscale)
    
        ax.bar(xvals[1], yvals[1], alpha=0.7, linewidth=0, color='black',
               align='center', clip_on=False, bottom=self.adjust)
        ax.bar(xvals[2], yvals[2], alpha=0.7, linewidth=0, color='orange',
               align='center', clip_on=False, bottom=self.adjust)
        ax.bar(xvals[3], yvals[3], alpha=0.7, linewidth=0, color='red',
               align='center', clip_on=False, bottom=self.adjust)
           
        # deal with the axis
        ax.axes.get_yaxis().set_visible(True)
        ax.tick_params(axis='y', direction='out', labelsize=6, left=True, right=False)
        ax.set_yticks( np.array(colthresh[1:])*heightscale+self.adjust )
        
        labels = [str(x) for x in colthresh[1:]]
        labels[-1] = '>'+labels[-1]
        ax.set_yticklabels( labels )
        ax.set_ylabel('Norm. React.', position=(0,0), ha='left', size=6)
           
        ax.set_frame_on(True)
        for l in ('right','top','bottom'):
            ax.spines[l].set_visible(False)
           
        ax.spines['left'].set_bounds(self.adjust, colthresh[3]*heightscale+self.adjust)


    
    def highlightLowReactivity(self, profile, cutoff=0.0001, window=1):
        """Profile is a np array (0-indexed) of values"""
        
        for i,r in enumerate(profile):
           
            if r<cutoff:
                patch = patches.Rectangle((i+0.5,0),window,4,linewidth=0,fill=True,alpha=0.2)
                self.topPatches.append(patch)     
                patch = patches.Rectangle((i+0.5,-4),window,4,linewidth=0,fill=True,alpha=0.2)
                self.botPatches.append(patch)     
                




    def writePlot(self, outPath="arcs.pdf", bounds=None, write=True,
                  msg=None, msg_pos=(0,1), msg_kwargs={}, **args):
        
        cutbounds = True
        if bounds is None:
            bounds = (0,len(self.seq))
            cutbounds = False
        else:
            bounds = [bounds[0]-0.5, bounds[1]+0.5] # make new copy

        
        doubleplot = len(self.botPatches)>0
        scaleFactor = 0.05
        
        figWidth = (bounds[1]-bounds[0])*scaleFactor
        figHeight = 2*max(self.height)*scaleFactor
        

        fig = plot.figure( figsize=(figWidth, figHeight)) # 500*scaleFactor
        fig.subplots_adjust(hspace=0.0)

        axT = fig.add_subplot(211)
        axB = None

        for pat in self.topPatches:
            axT.add_patch(pat)
            pat.set_clip_on(cutbounds)

        axT.set_xlim(bounds)
        axT.set_ylim(0, max(self.height))

        axT.set_frame_on(False)
        axT.axes.get_yaxis().set_visible(False)
        
        if self.toplegend:
            self.toplegend.drawLegend(axT, bounds[0], max(20, self.height[1]*.5))

        if doubleplot:
            axB = fig.add_subplot(212, sharex=axT)
            for pat in self.botPatches:
                axB.add_patch(pat)
                pat.set_clip_on(cutbounds)

            axB.set_xlim(bounds)
            axB.set_ylim(-max(self.height), 0)
            axB.set_frame_on(False)
            axB.axes.get_yaxis().set_visible(False)
            axB.tick_params(axis='both', which='both', top=False, bottom=False, labelbottom=False)
            
            if self.botlegend:
                self.botlegend.drawLegend(axB, bounds[0], -self.height[0]*.5) 



        # format sequence ticks
        axT.get_xaxis().tick_bottom()  
        
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(20)
        tickinterval = 5

        axT.xaxis.set_major_locator(majorLocator)
        axT.xaxis.set_major_formatter( FormatStrFormatter('%i') )
        axT.xaxis.set_minor_locator(minorLocator)
        axT.xaxis.set_minor_formatter( FormatStrFormatter('%i') )

        majlabels = axT.axes.get_xaxis().get_majorticklabels()
        minlabels = axT.axes.get_xaxis().get_minorticklabels()
        majticks = axT.axes.get_xaxis().get_major_ticks()
        minticks = axT.axes.get_xaxis().get_minor_ticks()
        

        for label in majlabels:
            label.set_size(10)
        
        if bounds[0] == 0: 
            majlabels[1].set_visible(False)
            minlabels[1].set_visible(False)
            majticks[1].set_visible(False)
            minticks[1].set_visible(False)
        
        beginint, r = divmod(bounds[0], 20)
        if r == 0:
            beginint+=1

        for i, label in enumerate(minlabels):
            label.set_size(7)
            if i % tickinterval == beginint:
                label.set_visible(False)


        # plot gridlines
        if self.grid:
            axT.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
            axT.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)

            if doubleplot:
                axB.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
                axB.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)



        # write title (structure names)
        if self.title:
            axT.text(0.0, self.height[1]-2, self.title, horizontalalignment='left',
                      size="24",weight="bold", color="blue")


        # put nuc sequence on axis 
        if self.drawseq:
            fontProp = matplotlib.font_manager.FontProperties(family = "sans-serif", 
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
            
            for i, nuc in enumerate(self.seq):
                
                if i<bounds[0]-1:
                    continue
                elif i>bounds[1]-1:
                    break

                if nuc == "T":
                    nuc = "U"

                try:
                    col = self.seqcolors[i]
                except IndexError:
                    col = "black"
                
                axT.annotate(nuc, xy=(i+0.5, 0),fontproperties=fontProp,color=col,
                             annotation_clip=False,verticalalignment="baseline")
        

        #if self.intdistance is not None:
        #    xvals = np.arange(bounds[0]+1, bounds[1]+1)
        #    yvals = self.intdistance[bounds[0]:bounds[1]+1]
        #    if np.mean(yvals) > 0:
        #        axT.plot( xvals, yvals, color="black", lw=2)
        #    else:
        #        axB.plot( xvals, yvals, color="black", lw=2)


        if self.reactprofile is not None:
            self.plotProfile(axT, bounds, **args)

        
        if msg is not None:
            axT.text(msg_pos[0], msg_pos[1], msg, transform=axT.transAxes, **msg_kwargs)


        # save the figure
        if write and outPath != 'show':
            fig.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)
            plot.close(fig)
        
        elif write and outPath == 'show':
            plot.show()

        else:
            return fig, axT, axB 



    def addCT(self, ctObj, color='black', alpha=0.7, panel=1):
    
        seqlen = len(ctObj.ct)

        if self.seq == '' or self.seq[0] == ' ':
            self.seq = ctObj.seq
        elif len(self.seq) != seqlen:
            print("Warning:: CT length = {0}; expecting length = {1}".format(seqlen, len(self.seq)))
        

        i = 0
        while i < seqlen:
            
            if ctObj.ct[i] > i:

                outerPair = [ i+1, ctObj.ct[i] ]

                # find the right side of helix
                lastPairedNuc = ctObj.ct[i]
                i += 1
                while i<seqlen and (lastPairedNuc-ctObj.ct[i]) == 1:
                    lastPairedNuc = ctObj.ct[i]
                    i += 1
                i -= 1


                innerPair = [ i+1, ctObj.ct[i] ]

                self.addArcPath(outerPair, innerPair, panel=panel, color=color, alpha=alpha)

            i += 1


    
    def addRings(self, ringfile, metric='z', panel=1, bins=None, contactfilter=(None,None),
                 filterneg=False):
        """Add arcs from ringmapper file
        metric  = z or sig
        """
        
        if contactfilter[0]:
            print("Performing contact filtering with distance={}".format(contactfilter[0]))
        
        if filterneg:
            print("Filtering out negative correlations")

        
        colors = [(44,123,182), (44,123,182), (171,217,233), (255,255,255), (253,174,97), (215,25,28), (215,25,28)]
        #colors = [(176,70,147), (176,70,147), (229,203,228), (255,255,255), (253,174,97), (215,25,28), (215,25,28)]



        colors = [tuple([y/255.0 for y in x]) for x in colors]
        
        alpha = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
        gradiant = [False, True, False, False, True, False]
        
        if metric=='z':
            if bins is None:
                bins = [-50, -5, -1, 0, 1, 5, 50]
            else:
                bins = [-1e10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]
        else:
            if bins is None:
                bins = [-1e10, -100, -20, 0, 20, 100, 1e10]
            else:
                bins = [-1e-10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]



        allcorrs = []
        
        with open(ringfile) as inp:
            header = inp.readline().split()
            window = int(header[1].split('=')[1])
            mname = header[2].split('=')[1]
            inp.readline()
                
            # define molecule length if not set from something else
            if self.seq == '':
                self.seq = ' '*int(header[0])
 

            for line in inp:
                spl = line.split()
                i = int(spl[0])
                j = int(spl[1])

                if contactfilter[0] and contactfilter[1].contactDistance(i,j) <= contactfilter[0]:
                    continue
                
                if filterneg and int(spl[3])<0:
                    continue

                if metric == 'z':
                    val = float(spl[4])
                else:
                    val = float(spl[2])

                if val > bins[3]:
                    allcorrs.append( ( int(spl[0]), int(spl[1]), float(spl[3])*val ) )


        if len(allcorrs)==0:
            print('WARNING: No RINGs passing filter in {}'.format(ringfile))


        for i in range(len(bins)-1):
            
            corrs = [c for c in allcorrs if bins[i]<c[2]<=bins[i+1]]
            
            corrs.sort(key=lambda x:x[2]) 
           
            for c in corrs:
                if gradiant[i]:
                    col = self._colorGrad(c[2], colors[i], colors[i+1], bins[i], bins[i+1])
                else:
                    col = colors[i]
                
                self.addArcPath( c[:2], panel=panel, color=col, alpha=alpha[i], window=window)

        
        # Add the legend
        t = 'RING {0} win={1} {2}'.format(mname, window, metric.upper())
        if filterneg:
            c = (colors[4], colors[5])
            l = ('>{}'.format(bins[4]), '>={}'.format(bins[5]))
        else:
            c = (colors[1], colors[2], colors[4], colors[5])
            l = ('<={}'.format(bins[1]), '<{}'.format(bins[2]), 
                 '>{}'.format(bins[4]), '>={}'.format(bins[5]))
        

        if panel>0:
            if self.toplegend is not None:
                self.toplegend.append(t, c, l)
            else:
                self.toplegend = ArcLegend(title=t, colors=c, labels=l)
        else:
            if self.botlegend is not None:
                self.botlegend.append(t, c, l)
            else:
                self.botlegend = ArcLegend(title=t, colors=c, labels=l)



    def addPairProb(self, dpObj, panel=1, bins=None):
        """ add pairing probability arcs from a dotplot object"""
 

        if self.seq == '':
            self.seq = ' '*dpObj.length
        elif len(self.seq) != dpObj.length:
            print("Warning:: dp file sequence length = {0}; CT/FASTA length = {1}".format(dpObj.length, len(self.seq)))


        refcolors = [ (150,150,150), (255,204,0),  (72,143,205) ,(81, 184, 72) ]
        refalpha = [0.55, 0.65, 0.75, 0.85]

        gradiant = [False, False, False, False]


        if bins is None:
            #bins = [0.01, 0.1, 0.5, 0.9, 1.0]
            bins = [0.03, 0.1, 0.3, 0.8, 1.0]
            colors = refcolors
            alpha = refalpha
        else:
            if len(bins) > 5:
                raise IndexError('{0} PairProb levels specified -- the maximum allowed is 4'.format(len(bins)-1))

            colors = refcolors[-(len(bins)-1):]
            alpha = refalpha[-(len(bins)-1):]


        # convert colors to percents
        colors = [tuple([y/255.0 for y in x]) for x in colors]
        

        for i in range(len(bins)-1):
            
            cutdp = dpObj.requireProb(bins[i], bins[i+1])
            cutpairs = cutdp.pairList()
            
            if len(cutpairs) == 0:
                continue
            
            # sort each of the pairs by its strength
            strength = []
            for p in cutpairs:
                elem = cutdp.partfun[p[0]-1]
                pos = ( elem['pair'] == p[1]-1) # make mask
                strength.append(elem['log10'][pos][0]) # get value

            strengthzip = list(zip(cutpairs, strength))
            strengthzip.sort(key=lambda x: x[1], reverse=True)

            cutpairs, strength = list(zip(*strengthzip))
            
            # iterate through sorted pairs and draw arcs
            for pindex, pair in enumerate(cutpairs):
 
                if gradiant[i]:
                    col = self._colorGrad(strength[pindex], colors[i], colors[i+1], 
                                          bins[i], bins[i+1], log=True)
                else:
                    col = colors[i]

                self.addArcPath(pair, panel=panel, color=col, alpha=alpha[i])


        # Add the legend
        t = 'Pairing Prob.'
        c = colors
        l = ['>{}'.format(x) for x in bins[:-1]]
        
        if panel>0:
            self.toplegend = ArcLegend(title=t, colors=c, labels=l)
        else:
            self.botlegend = ArcLegend(title=t, colors=c, labels=l) 
       




    def addPairMap(self, pmobj, panel=1, plotall=False):
        """plot pairs output by pairmapper
        plotall = True will plot all complementary correlations (not just 1&2)
        """

        if len(pmobj.primary)==0 and len(pmobj.secondary)==0:
            print('WARNING: No PAIR-MaP correlations in {}'.format(pmobj.sourcename))


        def getZ(corrlist):
            return [(x[0], x[1], x[3]) for x in corrlist]
        
        colors = [(100, 100, 100), (30,194,255), (0,0,243)]
        colors = [tuple([x/255.0 for x in c]) for c in colors]
 
        
        if len(self.seq)==0:
            self.seq = ' '*pmobj.molsize


        if plotall:
            self.plotAlphaGradient( getZ(pmobj.remainder), colors[0], (0.0,0.8),
                                    1, 6, window=pmobj.window, panel=panel)
        
        self.plotAlphaGradient( getZ(pmobj.secondary), colors[1], (0.2,0.6),
                                2, 6, window=pmobj.window, panel=panel)
        
        self.plotAlphaGradient( getZ(pmobj.primary), colors[2], (0.5,0.9),
                                2, 6, window=pmobj.window, panel=panel)
        
        
        # make legend
        c = colors[::-1]  # reverse colors so primary on top
        l = ['Principal','Minor','Not passing']
        if not plotall:
            c = c[:-1]
            l = l[:-1]
        
        if panel>0:
            if self.toplegend is not None:
                self.toplegend.append('PairMap', c, l)
            else:
                self.toplegend = ArcLegend(title='PairMap', colors=c, labels=l)
        else:
            if self.botlegend is not None:
                self.botlegend.append('PairMap', c, l)
            else:
                self.botlegend = ArcLegend(title='PairMap', colors=c, labels=l)



    def _colorGrad(self, value, colorMin, colorMax, minValue, maxValue, log=False):
        """ return an interpolated rgb color value """
        
        if log:
            value = 10**-value
            minValue = 10**-minValue
            maxValue = 10**-maxValue

        if value > maxValue:
            return colorMax
        elif value < minValue:
            return colorMin
        else:
            v = value - min(maxValue, minValue)
            v /= abs(maxValue-minValue)
        
            col = []
            for i in range(3):
                col.append( v*(colorMax[i] - colorMin[i]) + colorMin[i] )
            
            return col


    def plotAlphaGradient(self, pairlist, color, alpha_range, lb, ub, window=1, panel=1):
        """plot the list of pairs with the same color but different alpha"""

        if max(color) > 1:
            color = [x/255. for x in color]

        dynamicrange = ub-lb

        for i,j,val in sorted(pairlist, key=lambda x:x[2]):
            
            if val<lb:
                continue
            elif val>ub:
                scale = 1
            else:
                scale = (val-lb)/dynamicrange

            a = (alpha_range[1]-alpha_range[0])*scale+alpha_range[0]
            
            self.addArcPath( (i,j), window=window, color=color, panel=panel, alpha=a)



    def compareCTs(self, refCT, compCT, panel=1):
        
        if len(refCT.ct) != len(compCT.ct):
            raise IndexError('CT objects are different sizes')

        share = refCT.copy()
        refonly = refCT.copy()
        componly = refCT.copy()
        
        # ct files and set to blank (ie 0)
        for c in (share, refonly, componly):
            for i in range(len(c.ct)):
                c.ct[i] = 0
    

        for i in range(len(refCT.ct)):
            
            if refCT.ct[i] == compCT.ct[i]:
                share.ct[i] = refCT.ct[i]
            else:
                refonly.ct[i] = refCT.ct[i]
                componly.ct[i] = compCT.ct[i]
            
        sharedcolor = (150/255., 150/255., 150/255.)
        refcolor = (38/255., 202/255., 145/255.)
        compcolor = (153/255., 0.0, 1.0)
        self.addCT(share, color=sharedcolor, panel=panel)
        self.addCT(refonly, color=refcolor, panel=panel)
        self.addCT(componly, color=compcolor, panel=panel)
 

        sens,ppv,nums = refCT.computePPVSens(compCT, False)
        msg = 'Sens={0:.2f} PPV={1:.2f}'.format(sens, ppv)
        print(msg)

        if panel>0:
            self.toplegend = ArcLegend(colors=[sharedcolor,refcolor,compcolor], 
                                       labels=['Shared', refCT.name, compCT.name],
                                       msg=msg)                  
        else:
            self.botlegend = ArcLegend(colors=[sharedcolor,refcolor,compcolor], 
                                       labels=['Shared', refCT.name, compCT.name],
                                       msg=msg)                  
       
        return sens,ppv
 

    def assignSeqColors(self, shapearr):

        if len(shapearr) != len(self.seq):
            raise IndexError("Shape Array does not match length of sequence")

        col = []       
        for x in shapearr:

            if x < -4: 
                col.append( (160,160,160 ) )  # grey
            elif x > 0.85: 
                col.append( (255,0,0) )  # red
            elif 0.85 >= x >0.4: 
                col.append( (255,164,26) ) # orange
            else: 
                col.append( (0,0,0) ) # black

        self.seqcolors = [tuple([y/255.0 for y in x]) for x in col]


    
    def colorSeqByMAP(self, mapfile):
        
        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ')==len(self.seq):
            self.seq = seq
            
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.assignSeqColors(shape)
 
    
    def readSHAPE(self, mapfile):

        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ') == len(self.seq):
            self.seq = seq
        
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.reactprofile = shape
        
        # need to give the plot some height if no other things added
        if max(self.height)==0:
            self.height[1] = max(5, min(10, len(self.reactprofile)/50.))



    def readProfile(self, profilefile, dms=False):
        
        with open(profilefile) as inp:
            line = inp.readline().split()
            if len(line) == 2 or len(line)==4:
                ftype=1
            else:
                ftype=2
        

        if ftype==1:
            self.readSHAPE(profilefile)
        else:
            import ReactivityProfile
            profile = ReactivityProfile.ReactivityProfile(profilefile)
            self.reactprofile = profile.normprofile

            if self.seq == '' or self.seq.count(' ') == len(self.seq):
                self.seq = profile.sequence
        

        if dms:
            self.reactprofileType = 'DMS'
        
        # need to give the plot some height if no other things added
        if max(self.height)==0:
            self.height[1] = max(5, min(10, len(self.reactprofile)/50.))


        
    def old_addInteractionDistance(self, readMatFile, thresh, panel=1):
        
        height = []
        readMat = np.loadtxt(readMatFile)

        for i in range(len(self.seq)):
            
            for j in range(i,len(self.seq)):
                if readMat[i][j] < thresh:
                    break
            n = float(j-i)

            for j in range(i, -1, -1):
                if readMat[i][j] < thresh:
                    break
            m = float(i-j)

            try:
                z = n*math.sin(math.atan(m/n))
            except:
                z = 0

            height.append(panel*(z+self.adjust))

        self.intdistance = height



    def addInteractionDistance(self, readMatFile, thresh, panel=1):
        
        readmat = np.loadtxt(readMatFile)
        
        yvals = np.zeros(len(self.seq))

        for i in range(len(self.seq)):

            for j in range(i, len(self.seq)):
                if readmat[i,j] >=thresh:

                    idx = int(round((j-i+1)/2.0))
                    
                    if idx > yvals[i+idx]:
                        yvals[i+idx] = idx
        
        self.intdistance = panel*yvals


#############################################################################




def parseArgs():

    prs = argparse.ArgumentParser()
    prs.add_argument("outputPDF",type=str, help="Name of the output PDF graphic")
    prs.add_argument("--ct", type=str, help="Base CT file to plot. By default, the first (0th) structure is plotted. Alternative structures can be selected by passing ,# after the CT file name (e.g. --ct ctfile.ct,3)")
    prs.add_argument("--fasta", type=str, help="Fasta sequence file")
    prs.add_argument("--refct", type=str, help="Reference CT to compare base ct against. By default the first (0th) structure is plotted. Alternative structures can be selected by passing ,# after the CT file name (e.g. --refct ctfile.ct,3)")
    prs.add_argument("--probability", type=str, help="Pairing probability file in dotplot format. By default, arcs are drawn for the following probability intervals: [0.03,0.1], [0.1,0.3], [0.3,0.8], [0.8,1.0]. These intervals can be modified by passing thresholds as comma-separated values. For example, --prob file.dp,0.03,0.1,0.3,0.8,1.0 specifies the default intervals. At most 4 intervals can be plotted, but fewer intervals are allowed (e.g. --prob file.dp,0.1,0.3,0.8,1.0 will plot 3 intervals).")
    
    prs.add_argument("--ringz", type=str, help="Plot Z-scores from ringmapper correlation file. Default color thresholds are Z=1 and Z=5. Can be modified by passing thresholds as comma-separated values (e.g. --corrz corrs.txt,1,5)")
    
    prs.add_argument("--ringsig", type=str, help="Plot statistical significance from ringmapper file. Default color thresholds are [20,500]. Can be modified by passing thresholds as comma-separated values (e.g. --corrsig corrs.txt,20,500)")

    prs.add_argument("--pairmap", type=str, help="Plot pairmap signals from pairmap file. By default plot principal & minor correlations. Can plot all complementary correlations by passing ,all (e.g. --pairmap pairmap.txt,all)")

    prs.add_argument("--compare_pairmap", type=str, help="Plot pairmap signals from second pairmap file. By default plot principal & minor correlations. Can plot all complementary correlations by passing ,all (e.g. --compare_pairmap pairmap.txt,all)")


    prs.add_argument("--ntshape", type=str, help="Color nucs by shape reactivty in provided shape/map file")  
    prs.add_argument("--profile", type=str, help="Plot reactivity profile on top from shape/map file")
    
    prs.add_argument("--dmsprofile",type=str, help='Normalize and plot DMS reactivity from profile file')

    prs.add_argument("--bottom", action='store_true', help="Plot arcs on the bottom")

    prs.add_argument("--title", type=str, default='', help='Figure title')
    prs.add_argument("--showGrid", action="store_true", help="plot a grid on axis ticks")
    #prs.add_argument("--intDistance", type=str, help="Pass depth file for computing likely max interaction distance")
    #prs.add_argument("--depthThresh", type=int,default=10000, help="Interaction depth threshold (Default=10,000)")
    prs.add_argument("--bound", type=str, help="comma separated bounds of region to plot (e.g. --bound 511,796)")

    prs.add_argument("--filternc", action="store_true", help="filter out non-canonical pairs in ct")
    prs.add_argument("--filtersingle", action="store_true", help="filter out singleton pairs in ct")

    prs.add_argument("--contactfilter", type=int, help="filter rings by specified contact distance (int value)")
    
    prs.add_argument("--filternegcorrs", action="store_true", help='filter out negative correlations')  

    args = prs.parse_args()
 

    if args.refct and not args.ct:
        exit("--refct is invalid without --ct")


    numplots = int(args.ct is not None)


    # subparse the prob argument
    args.probability_bins = None
    if args.probability:
        numplots += 1

        spl = args.probability.split(',')
        try:
            if len(spl) > 1:
                args.probability_bins = [float(x) for x in spl[1:]]
                args.probability = spl[0]
        except:
            raise TypeError('Incorrectly formatted --probability argument {}'.format(args.probability))
    
    

    def subparse3(arg, name):
        
        outarg = arg
        bins = None

        spl = arg.split(',')
        try:
            if len(spl) == 3:
                bins = [float(x) for x in spl[1:]]
                outarg = spl[0]
            elif len(spl) != 1:
                raise TypeError
        except:
            raise TypeError('Incorrectly formatted {0} argument {1}'.format(name, arg))

        return outarg, bins


    if args.ringz:
        numplots += 1
        args.ringz, args.ringz_bins = subparse3(args.ringz, '--ringz')
    
    if args.ringsig:
        numplots += 1
        args.ringsig, args.ringsig_bins = subparse3(args.ringsig, '--ringsig')


    args.pairmap_all = False
    if args.pairmap:
        numplots += 1
        spl = args.pairmap.split(',')
        if len(spl)==1:
            pass
        elif len(spl)==2 and spl[1] == 'all':
            args.pairmap_all = True
            args.pairmap = spl[0]
        else:
            raise TypeError('Incorrectly formatted --pairmap argument {}'.format(args.pairmap))
    

    args.ringpairsuper = False
    if args.pairmap and (args.ringz or args.ringsig) and args.ct:
        args.ringpairsuper = True

    if numplots > 2 and not args.ringpairsuper:
        exit('Too many plots! Please select at maximum 2 of [--ct, --probability, --ringz, --ringsig, --pairmap --compare_pairmap]')


    if args.contactfilter and not args.ct:
        exit('Cannot perform contact filtering without --ct file')

    if args.ringpairsuper and args.contactfilter is None:
        args.contactfilter = 20

    # subparse the bounds argument
    if args.bound:
        args.bound = map(int, args.bound.split(','))
 
    
    # subparse the ct arguments
    if args.ct:
        spl = args.ct.split(',')
        if len(spl)==1:
            args.ctstructnum = 0
        else:
            args.ctstructnum = int(spl[1])
            args.ct = spl[0]

    if args.refct:
        spl = args.refct.split(',')
        if len(spl)==1:
            args.refctstructnum = 0
        else:
            args.refctstructnum = int(spl[1])
            args.refct = spl[0]


    return args





if __name__=="__main__":
 
    args = parseArgs()
    
    msg = None
    CT1=None
    
    aplot = ArcPlot(title = args.title, fasta=args.fasta)

    panel = 1
    if args.bottom:
        panel = -1

    if args.ct:
        
        CT1 = RNAtools.CT(args.ct, structNum=args.ctstructnum, filterNC=args.filternc, filterSingle=args.filtersingle)

        if args.refct:
            refCT = RNAtools.CT(args.refct, structNum=args.refctstructnum, filterNC=args.filternc, filterSingle=args.filtersingle)
            aplot.compareCTs( refCT, CT1, panel=panel)
        
        else:
            alpha = 0.7
            if args.ringpairsuper:
                alpha=0.2
            
            aplot.addCT( CT1, panel=panel, alpha=alpha)

        panel *= -1


    if args.probability:
        aplot.addPairProb( RNAtools.DotPlot(args.probability), panel=panel, bins=args.probability_bins)
        panel *= -1


    if args.pairmap:

        from pmanalysis import PairMap
        
        if args.ringpairsuper:
            aplot.addPairMap( PairMap(args.pairmap), panel=1, plotall=args.pairmap_all)

        else:
            aplot.addPairMap( PairMap(args.pairmap), panel=panel, plotall=args.pairmap_all)
            panel *= -1
        

    if args.ringz:
        aplot.addRings(args.ringz, panel=panel, metric='z', bins=args.ringz_bins,
                       contactfilter=(args.contactfilter, CT1), filterneg=args.filternegcorrs)
        panel *= -1
    
    if args.ringsig:
        aplot.addRings(args.ringsig, panel=panel, metric='sig', bins=args.ringsig_bins,
                       contactfilter=(args.contactfilter, CT1), filterneg=args.filternegcorrs)
        panel *= -1

    
    if args.compare_pairmap:
        
        from pmanalysis import PairMap
        
        aplot.addPairMap( PairMap(args.compare_pairmap), panel=panel, plotall=args.pairmap_all)
        panel *= -1




    #if arg.intDistance:
    #    aplot.addInteractionDistance(arg.intDistance, arg.depthThresh, panel)
            

    if args.ntshape:
        aplot.colorSeqByMAP(args.ntshape)

    if args.profile:  
        aplot.readProfile(args.profile)
    
    if args.dmsprofile:
        aplot.readProfile(args.dmsprofile, dms=True)

    if args.showGrid:
        aplot.grid = True

    aplot.writePlot( args.outputPDF, bounds = args.bound, msg=msg)







# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
