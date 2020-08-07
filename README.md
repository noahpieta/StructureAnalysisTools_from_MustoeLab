# RNATools
Houses the RNATools script started by Gregg and now altered by like everybody. 

# foldPK
Functions and scripts that wrap around ShapeKnots to fold RNAs with multiple PKs.
Structure models are written as CT format files (see RNAstructure documentation for details). 
Multiple CT files may be generated as part of the hierarchical folding process. 
These are denoted as [outprefix].1.ct, .2.ct, etc. The final solution will be named [outprefix].f.ct.
A constraint file may also be written as [outprefix].cons during folding; this can be deleted.

Note that if calling foldPK from the command line UI you will either need to use the --skpath
argument to pass the path to the ShapeKnots executable OR edit the 'default' path located 
within the code (line 103) to point to your ShapeKnots executable
