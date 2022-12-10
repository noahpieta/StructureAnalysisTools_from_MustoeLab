Tools for analyzing and visualizing RNA probing data and structure models

Copyright 2022 Anthony Mustoe

Contact: anthony.mustoe@bcm.edu
-------------------------------


# ArcPlot

Description
-----------
Can plot the following types of data

-Structures
    --ct (plot a single ct file)
    --refct (Compare ct to a reference ct and plot both)

-Pairing probability
    --pairprob (plot pairing probability from dotplot file)

-Rings
    --ringz (plot z-scores)
    --ringsig (plot correlation Chi2 significance)

-PAIR-MaP
    --pairmap (plot PAIR-MaP base pairs)

-Reactivity profiles
    --profile (plot normalized SHAPE reactivity data)
    --dmsprofile (plot normalized reactivity data colored to slightly different scale)


# StructureObjects
Houses CT and DotPlot objects. Renamed from RNATools


# foldPK
Functions and scripts that wrap around ShapeKnots to fold RNAs with multiple PKs.
Structure models are written as CT format files (see RNAstructure documentation for details). 
Multiple CT files may be generated as part of the hierarchical folding process. 
These are denoted as [outprefix].1.ct, .2.ct, etc. The final solution will be named [outprefix].f.ct.
A constraint file may also be written as [outprefix].cons during folding; this can be deleted.

Note that if calling foldPK from the command line UI you will either need to use the --skpath
argument to pass the path to the ShapeKnots executable OR edit the 'default' path located 
within the code (line 103) to point to your ShapeKnots executable


