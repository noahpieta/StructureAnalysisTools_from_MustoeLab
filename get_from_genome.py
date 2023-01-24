import pandas as pd
import pybedtools
#from pyliftover import LiftOver
import argparse
import os

class APtools():
    
    def __init__(self):
        self.start_coors = []
        self.end_coors = []
        
    def get_rna_info(self, bedfile):
        self.exons = [0]
        
        with open(bedfile, 'r') as f:
            for line in f:
                spl = line.rstrip().split()
                self.chr = spl[0]
                self.start_coors.append(int(spl[1]))
                self.end_coors.append(int(spl[2]))
                self.strand = spl[3]
                
                exon = int(spl[2]) - int(spl[1])
                self.exons.append(exon) # use this exon length to update transcript coordinate
         
        self.exons = self.exons[:-1]
        
        self.start = min(self.start_coors)
        self.end = min(self.end_coors)
        self.length = self.end - self.start + 1
        
        if 'chr' not in self.chr:
            self.chr = 'chr' + self.chr
            
    def convert_coordinate(self, genome1, genome2):
        lo = LiftOver(genome1, genome2)
        self.start_coors2 = [] 
        self.end_coors2 = []
        
        for _start, _end in zip(self.start_coors, self.end_coors):
            start = lo.convert_coordinate(self.chr, _start, strand = self.strand)[0][1]
            self.start_coors2.append(start)

            end = lo.convert_coordinate(self.chr, _end, strand = self.strand)[0][1]
            self.end_coors2.append(end)

    # create fasta file from genomic coordinate using pybedtools get fasta function
    def get_fasta(self, chrom, start_coors, end_coors, strand, genome, fastaname):

        seq = ''
        for start, end in zip(self.start_coors, self.end_coors):
            chrom = self.chr.replace('chr', '')
            line = f'{chrom} {start} {end} . . {strand}'
            bt = pybedtools.BedTool(line, from_string = True)
            try:
                _seq = bt.sequence(fi = genome, s = True)
                _seq = open(_seq.seqfn).readlines()[1].rstrip()
                seq += _seq
            except:
                print('Provided coordinate(s) does not match with a sequence in the reference genome file. Exiting.')
                exit()

        name = fastaname.split('.')[0]
        with open(fastaname,'w') as fa:
            fa.write(f'>{name}\n{seq.upper()}')

        print(f'Sequence is in {fastaname}')
        return fastaname
    
    # create reactivity file from genomic coordinate 
    def get_reactivity(self, chrom, start_coors, end_coors, strand, reactivity_file, outname):
        if reactivity_file.endswith('gz'):
            compression = 'gzip'
        else:
            compression = None
            
        reactivity = pd.read_csv(reactivity_file, sep = '\t', compression = compression, \
                         names = ['chr', 'start', 'end', 'NA', 'score', 'strand'])
        
        filtered = pd.DataFrame()
        length = 0
        scores = []
        
        #lo = LiftOver('hg19', 'hg38')
        
        for start, end in zip(start_coors, end_coors):
            temp = reactivity[(reactivity['chr']==chrom)&(reactivity['strand']==strand)]
            temp = temp[(start<=temp['start'])&(temp['start']<=end)]
            filtered = pd.concat([filtered, temp])
            _length = end - start + 1
            length += _length
            
            for i in range(start, end+1):
                try:
                    score = filtered.loc[filtered['start']==i, 'score'].values[0]
                except:
                    score = -999
            
                scores.append(score)
                
        if strand == '-':
            scores.reverse()
            
        with open(outname, 'w') as out:
            for idx, score in enumerate(scores, start = 1):
                out.write(f'{idx}\t{score}\n')
                
        print(f'Reactivity is in {outname}')
        
        return outname

    def convert_ucsc(self, ucsc, chrom, start_coors, end_coors, strand):

        with open(ucsc, 'r') as f:
            for idx, line in enumerate(f):
               if line.startswith('variable') or line.startswith('fixed'):
                   skip = idx + 1
                   break

        df = pd.read_csv(ucsc, sep='\t', skiprows = skip, names = ['pos', 'score'])
        start = int(df.iloc[0]['pos'])
        end = int(df.iloc[-1]['pos'])
        
        filtered = pd.DataFrame()
        length = 0
        scores = []
        
        #lo = LiftOver('hg19', 'hg38')
        
        for start, end in zip(start_coors, end_coors):
            temp = df[(start<=df['pos'])&(df['pos']<=end)]
            filtered = pd.concat([filtered, temp])
            _length = end - start + 1
            length += _length
           
            for i in range(start, end+1):
                try:
                    score = filtered.loc[filtered['pos']==i, 'score'].values[0]
                except:
                    score = score
            
                scores.append(score)

        # reverse list of scores if antisense strand	
        if strand == '-':
            scores.reverse()
        
        scores = scores[:-1]
        name = os.path.basename(ucsc)
        basename, ext = os.path.splitext(name)

        with open(f'{basename}_parsed.txt', 'w') as out:
            for idx, score in enumerate(scores, start = 1):
                out.write(f'{idx}\t{score}\n')
                
        print(f'Conservation is in {basename}_parsed.txt')
 
#########################################################################################

def parseArgs():
    prs = argparse.ArgumentParser()
    prs.add_argument("--coordinate", type=str, help='Genomic coordinate in BED format for RNA of interest. \
                                                Require 4 columns: chromosome start end strand. \
                                                Each line corresponds to one exon.')
    prs.add_argument("--genome", type=str, help='Path to genome assembly to retrieve sequence from.')
    prs.add_argument("--fasta_out", type=str, help='Name of output FASTA file.')
    prs.add_argument("--reactivity", type=str, help='Path to transcriptome wide reactivtity in BED format.')
    prs.add_argument("--reactivity_out", type=str, help='Name of output reactivity file (SHAPE or DMS).')
    prs.add_argument("--liftFA", type=str, help='Use USCS LiftOver utilities to convert coordinates before \
                                                generating a FASTA file. Input 2 genomes separated by comma.')
    prs.add_argument("--liftRE", type=str, help='Use USCS LiftOver utilities to convert coordinates before \
                                                generating a reactivity file. Input 2 genomes separated by comma.')
    prs.add_argument("--ucsctrack",type=str, help='File containing UCSC track.')
    args = prs.parse_args()
    return args

#########################################################################################

if __name__=="__main__":
 
    args = parseArgs()
    
    APt = APtools()
    
    if args.coordinate:
        APt.get_rna_info(args.coordinate)
        
    # get fasta arguments
    if args.genome:
        if args.liftFA:
            genomes = args.liftFA.split(',')
            APt.convert_coordinate(genomes[0], genomes[1])
            APt.get_fasta(APt.chr, APt.start_coors2, APt.end_coors2, APt.strand, args.genome, args.fasta_out)
        else:
            APt.get_fasta(APt.chr, APt.start_coors, APt.end_coors, APt.strand, args.genome, args.fasta_out)
            
    # get reactivity arguments
    if args.reactivity:
        if args.liftRE:
            genomes = args.liftRE.split(',')
            APt.convert_coordinate(genomes[0], genomes[1])
            APt.get_reactivity(APt.chr, APt.start_coors2, APt.end_coors2, APt.strand, args.reactivity, args.reactivity_out)
        else:
            APt.get_reactivity(APt.chr, APt.start_coors, APt.end_coors, APt.strand, args.reactivity, args.reactivity_out)

    # get conservation arguments
    if args.ucsctrack:
        APt.convert_ucsc(args.ucsctrack, APt.chr, APt.start_coors, APt.end_coors, APt.strand)
