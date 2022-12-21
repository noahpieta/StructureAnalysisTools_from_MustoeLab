# Anthony Mustoe
# 2020

import RNAStructureObjects as RNAtools
import subprocess, argparse, os



def writeConstraints(outfile, pairs):
    """pairs is a list of base paired nucs to restrain to single-stranded
    for further pk folding""" 
   
    # turn pairs into 1D list of nts
    nts = [i for p in pairs for i in p]
    nts.sort()

    with open(outfile, 'w') as out:

        out.write('DS:\n-1\nSS:\n')
        
        for n in nts:
            out.write("{}\n".format(n))
        out.write('-1\nMod:\n-1\nPairs:\n-1 -1\nFMN:\n-1\nForbids:\n-1 -1\n')

 

def runShapeKnots(ShapeKnotsPath, seqfile, output, shapefile=None, shapepars=None, dmsfile=None, bpfile=None, consfile=None, cordero=False):

    commands = [ShapeKnotsPath, '-m','1']

    if shapefile is not None and dmsfile is not None:
        raise AttributeError('Only dms or shape can be specified')

    elif shapefile is not None:
        commands.extend(['-sh', shapefile])
    
    if shapepars is not None:
        commands.extend(['-sm',str(shapepars[0]), '-si',str(shapepars[1])])


    if dmsfile is not None and not cordero:
        commands.extend(['-dmsnt',dmsfile])
    elif dmsfile is not None and cordero:
        commands.extend(['-dms',dmsfile])
    
    if bpfile is not None:
        commands.extend(['-x', bpfile])
    if consfile is not None:
        commands.extend(['-c', consfile])
    
    commands.extend([seqfile, output])

    print(commands)
    subprocess.call(commands)




def iterativeShapeKnots(ShapeKnotsPath, seqfile, outprefix, shapefile=None, shapepars=None, dmsfile=None, bpfile=None, maxPKs=5, cordero=False):

    modelnumber = 0
    hasPK = True
    pkpairs = []
    
    consfile = None
    consfilename = '{}.cons'.format(outprefix)

    while hasPK and modelnumber < maxPKs:
        
        modelnumber += 1
        
        if len(pkpairs) == 0:
            consfile = None
        else:
            writeConstraints(consfilename, pkpairs)
            consfile = consfilename

        output = '{0}.{1}.ct'.format(outprefix, modelnumber)
        runShapeKnots(ShapeKnotsPath, seqfile, output, shapefile=shapefile, shapepars=shapepars, dmsfile=dmsfile, bpfile=bpfile, 
                      consfile=consfile)

        ct = RNAtools.CT(output)
        pk = ct.extractPK()
        
        if pk[0] is None:
            hasPK = False
        else: 
            print("PK found! {}".format(pk[0]))
            pkpairs.extend(pk[0])  
        

    # add pairs back in
    if len(pkpairs)>0:
        ct.addPairs(pkpairs)
        ct.writeCT('{}.f.ct'.format(outprefix))
    else:
        subprocess.call(['mv', output, '{}.f.ct'.format(outprefix)])





if __name__ == '__main__':
    
    prs = argparse.ArgumentParser(description='Script to iteratively run ShapeKnots to find multiple PKs')
    prs.add_argument('seqfile', help='molecule name')
    prs.add_argument('outprefix', help="prefix to append to ct files. Note that .f.ct is the 'final' fold")
    prs.add_argument('--dmsfile', help='dms reactivity file')
    prs.add_argument('--shapefile', help='shape reactivity file')
    prs.add_argument('--sm', help='shape slope (sm) parameter if different than ShapeKnots default')
    prs.add_argument('--si', help='shape intercept (si) parameter if different than ShapeKnots default')
    prs.add_argument('--bpfile', help='Pairmapper bp restraint file')
    prs.add_argument('--skpath', help='Path to ShapeKnots executable. If not defined, will try and use internally defined path. You can also modify this internal path if desired (will have to edit within code)')
    args = prs.parse_args()
    
    
    if (args.sm is None) + (args.si is None) == 1:
        exit('Need to specify both --sm and --si')

    if args.sm is not None and args.si is not None:
        args.shapepars = (float(args.sm), float(args.si))
    else:
        args.shapepars = None


    # default path to ShapeKnots. Change to path on your Machine
    defaultpath = '/Users/anthonymustoe/Code/RNAstructure/exe/ShapeKnots'
     

    if args.skpath is not None and not os.path.isfile(args.skpath):
        exit('User-provided path to ShapeKnots is invalid! File does not exist')
    elif args.skpath is None:
        args.skpath = defaultpath
        if not os.path.isfile(args.skpath):
            exit('Default path to ShapeKnots {} is invalid! File does not exist'.format(args.skpath))


    iterativeShapeKnots(args.skpath, args.seqfile, args.outprefix, 
                        shapefile=args.shapefile, shapepars=args.shapepars,
                        dmsfile=args.dmsfile, bpfile=args.bpfile)
    
    
