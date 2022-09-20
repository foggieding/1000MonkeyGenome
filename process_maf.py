## Script: get_invMaf_aln.py
## Description: infer inversion status from MAF alignment file.
## Author: Kevin Lee
## Date: 2022.01.06

# NOTE:
# (1) Inversion name must be formated as starting with chromosome and delimited with "_", e.g, chr19_inv1 or chr19_SV_12

import argparse
from utils import InputStream
from utils import PathType
import sys,os
import numpy as np
import re
import pickle

def warn(message):
    prog = os.path.basename(sys.argv[0])
    sys.stderr.write(prog + ": " + message + "\n")

def mafBlocks(beg1, beg2, seq1, seq2):
    '''Get the gapless blocks of an alignment, from MAF format.'''
    '''mafBlocks()输出的内容为maf中没有gap的1to1 alignment的interval，格式是物种1的start位点，物种2的start位点以及物种1和物种2无gap alignment的size'''
    size = 0
    for x, y in zip(seq1, seq2):
        if x == "-":
            if size:
                yield beg1, beg2, size
                beg1 += size
                beg2 += size
                size = 0
            beg2 += 1
        elif y == "-":
            if size:
                yield beg1, beg2, size
                beg1 += size
                beg2 += size
                size = 0
            beg1 += 1
        else:
            size += 1
    if size: yield beg1, beg2, size

def alignmentInput(stream):
    '''Get alignments and sequence lengths, from MAF or tabular format.'''
    mafCount = 0
    for line in stream:
        #print("I'm line in maf.gz: ", line)
        w = line.split()
        if not w:continue

        #print("I'm w[0],see whether I'm s: ",w[0])
        #print("I'm line[0], see How I look like: ",line[0])

        if w[0] == 's':
            #print("I'm starting with s")
            if mafCount == 0:
                 #print("I'm in mafCount 0 ")
                chr1, beg1, seqlen1, strand1, seq1 = w[1], int(w[2]), int(w[5]), w[4], w[6]
                if w[4] == "-": beg1 -= seqlen1
                mafCount = 1
            else:
                #print("I'm in mafCount 1 ")
                chr2, beg2, seqlen2, strand2,seq2 = w[1], int(w[2]), int(w[5]), w[4], w[6]
                if w[4] == "-": beg2 -= seqlen2
                blocks = list(mafBlocks(beg1, beg2, seq1, seq2)) # 生成物种1和物种2每一条maf alignment之间的一一对应的block
                #print(chr1, strand1, chr2, strand2, blocks)
                yield chr1, strand1, chr2, strand2, blocks  # 返回物种1染色体（hg38.chr1），物种1比对的方向，物种1染色体（calJac4.chr9），物种2比对的方向，物种1和物种2一一对应的mafblock

                mafCount = 0

def parseChrSize(args):
    with InputStream(args.chrSize) as stream:
        for chr in stream:
            chrl = chr.rstrip().split('\t')
            yield chrl[0],int(chrl[1])

def flankRanges(chr,spos,epos,chrSize,args):
    size = epos - spos
    spos_f1 = max(0,spos - size * args.flankFrac)
    chr = re.sub('^.*\\.','',chr)
    epos_f2 = min(chrSize[chr],epos + size * args.flankFrac)
    return spos_f1,spos,epos,epos_f2

def readInv(args,chrSize):
    # Read inversion list
    # Note: generate flanking region as well
    with InputStream(args.inv) as stream:
        for inv in stream:
            # inversion itself
            invl = inv.rstrip().split('\t')
            #yield invl[0],invl[1],invl[2],invl[3]
            # flanking intervals of inversion
            spos_f1,epos_f1,spos_f2,epos_f2 = flankRanges(invl[0],int(invl[1]),int(invl[2]),chrSize,args)
            yield invl[3],((int(invl[1]),int(invl[2])),(spos_f1,epos_f1),(spos_f2,epos_f2))

def processMaf(args):
    with InputStream(args.maf) as stream:
        #print(stream)
        mafRecords = list(alignmentInput(stream)) # 格式：[(hg38.chr1, "+", calJac4.chr9, "-", [(hg38_start,calJac4_start,size),(hg38_start,calJac4_start,size)]),(...)]
    #print("mafRecords: ",mafRecords)
    with open(args.processMaf,'wb') as f:
        pickle.dump(mafRecords,f)

    return mafRecords

def determineInvOri():
    warn('determining inversion orientation')

def cropMaf(invRange,blocks):
    # 请确保物种1中的maf alignment方向全部为+。因为物种1 maf alignment方向为-（即headBeg1为-的情况未在考虑范围内）
    #warn('entering cropMaf()')
    headBeg1, headBeg2, headSize = blocks[0]
    invSpos,invEpos = invRange
    if headBeg1 < 0:
        invSpos, invEpos = -invEpos, -invSpos
    for beg1, beg2, size in blocks:
        b1 = max(invSpos, beg1)
        e1 = min(invEpos, beg1 + size)
        if b1 >= e1: continue
        offset = beg2 - beg1
        b2 = b1 + offset
        e2 = e1 + offset
        if b2 >= e2: continue
        yield b1, b2, e2 - b2

def intersectInvMaf(invl,mafRecords,args):
    warn('going into intersectInvMaf')
    intersect1 = dict() # upstream region
    intersect2 = dict() # downstream region
    intersect = dict() # inv itself

    for inv in invl:
        intersect1[inv[0]] = []
        intersect2[inv[0]] = []
        intersect[inv[0]] = []
        for maf in mafRecords:
            ## check if maf records locate in the same chromosome as inv
            invChr = inv[0].split('_')[0]
            invChr = invChr[0].lower() + invChr[1:] # lowercase first character to be compatible with strand-seq inversion name
            mafChr1 = re.sub('.*?\.','',maf[0])
            if invChr != mafChr1:
                continue

            ## get chr of inv and maf
            invName,invSpos,invEpos,invSpos1,invEpos1,invSpos2,invEpos2 = inv[0],inv[1][0][0],inv[1][0][1],inv[1][1][0],inv[1][1][1],inv[1][2][0],inv[1][2][1]
            mafStrand1,mafChr2,mafStrand2,blocks  = maf[1],re.sub('.*?\.','',maf[2]),maf[3],maf[4]

            ## crop maf given inv
            # For region upstream inv, inv itself, region downstream inv
            cropRegion1 = list(cropMaf((invSpos1,invEpos1),blocks))
            cropRegion2 = list(cropMaf((invSpos2,invEpos2),blocks))
            cropRegion = list(cropMaf((invSpos,invEpos),blocks))
            if cropRegion1:
                intersect1[invName].append(tuple([mafChr1,mafStrand1,mafChr2,mafStrand2,list(cropMaf((invSpos1,invEpos1),blocks))]))
            if cropRegion2:
                intersect2[invName].append(tuple([mafChr1,mafStrand1,mafChr2,mafStrand2,list(cropMaf((invSpos2,invEpos2),blocks))]))
            if cropRegion:
                intersect[invName].append(tuple([mafChr1,mafStrand1,mafChr2,mafStrand2,list(cropMaf((invSpos,invEpos),blocks))]))

    print('cropRegion: ',intersect,file=sys.stderr,flush=True)
    ## save intersect list
    with open(args.intersect,'wb') as f:
        pickle.dump(intersect,f)
    with open(args.intersect1,'wb') as f:
        pickle.dump(intersect1,f)
    with open(args.intersect2,'wb') as f:
        pickle.dump(intersect2,f)

def runFromArgs(args):
    if args.inv:
        print('Handling inversion list...',file=sys.stderr,flush=True)
        # parse chromosome size
        chrSize = dict(parseChrSize(args))
        #print(chrSize)
        # handle inversion list
        invl = list(readInv(args,chrSize)) # 格式: flank_f1, inversion_itself, flank_f2 [('chr1_inv1', ((1007565, 1009112), (1004471.0, 1007565), (1009112, 1012206.0))), ('chr1_inv2', ((1383839, 1384532), (1382453.0, 1383839), (1384532, 1385918.0)))]
        print(invl)
        if args.invlst:
            with open(args.invlst,'wb') as f:
                pickle.dump(invl,f)
    elif args.invlst:
        print('Read invlst',file=sys.stderr,flush=True)
        with open(args.invlst, "rb") as f:  # Unpickling
            invl = pickle.load(f)
    else:
        print('Please specify inversion list to run from scratch or provide pickled inversions',file=sys.stderr,flush=True)
        sys.exit(0)

    if args.invOnly:
        print('Only pickle inversion list...Exit...',file=sys.stderr,flush=True)
        sys.exit(0)
    else:
        print('Processing MAF alignment...',file=sys.stderr,flush=True)
        #with open(args.test,'w') as f:
        #    f.write('a\n')
        mafRecords = processMaf(args) # 格式：[(hg38.chr1, "+", calJac4.chr9, "-", [(hg38_start,calJac4_start,size),(hg38_start,calJac4_start,size)]),(...)]
        #print(mafRecords)

        print('Intersect inversion and flanking regions with MAF alignment',file=sys.stderr,flush=True)
        intersectInvMaf(invl,mafRecords,args)

    '''
    if args.inv:
        warn('Handling inversion list...')
        # parse chromosome size
        chrSize = dict(parseChrSize(args))
        #print(chrSize)
        # handle inversion list
        invl = list(readInv(args,chrSize)) # 格式: flank_f1, inversion_itself, flank_f2 [('chr1_inv1', ((1007565, 1009112), (1004471.0, 1007565), (1009112, 1012206.0))), ('chr1_inv2', ((1383839, 1384532), (1382453.0, 1383839), (1384532, 1385918.0)))]
        print(invl)
        if args.invlst:
            with open(args.invlst,'wb') as f:
                pickle.dump(invl,f)
    elif args.invlst:
        warn('Read invlst')
        with open(args.invlst, "rb") as f:  # Unpickling
            invl = pickle.load(f)
    else:
        warn('Please specify inversion list to run from scratch or provide pickled inversions')
        sys.exit(0)

    if args.invOnly:
        warn('Only pickle inversion list...Exit...')
        sys.exit(0)
    else:
        warn('Processing MAF alignment...')
        mafRecords = processMaf(args) # 格式：[(hg38.chr1, "+", calJac4.chr9, "-", [(hg38_start,calJac4_start,size),(hg38_start,calJac4_start,size)]),(...)]
        #print(mafRecords)

        warn('Intersect inversion and flanking regions with MAF alignment')
        intersectInvMaf(invl,mafRecords,args)
    '''

if __name__ == "__main__":
    descrition='Infer inversion status from MAF alignment file.'
    epilog='The input MAF file may be gzipped. If the input file is omitted then input is read from stdin. Output is written to stdout.'
    parser = argparse.ArgumentParser(description=descrition,epilog=epilog)
    parser.add_argument('--inv',metavar='<BED>',default=None,help='Inversion list in bed file format')
    parser.add_argument('--invOnly',default=False,action='store_true',help="If specified, only pickle inversion list")
    parser.add_argument('--test',default=None,help='test')
    parser.add_argument('--maf',metavar='<GEM>',default=None,help='MAF alignment')
    parser.add_argument('--flankFrac',metavar='<FRAC>',default=2,type=np.float32,help='flank length used to determine background orientation')
    parser.add_argument('--chrSize',metavar='<TXT>',default=None,help='chromosome size file for species1')
    parser.add_argument('--outDir',metavar='<DIR>',default='.', type=PathType(exists=True, type='dir'), help='Output directory')
    parser.add_argument('--invlst', metavar='<PDS>', default=None, help='File to save inversion-related intervals and will be served as input when --inv not assigned.')
    parser.add_argument('--intersect1',metavar='<PDS>', help='File to save intersect between upsteram region of inv and maf alignment')
    parser.add_argument('--intersect2',metavar='<PDS>', help='File to save intersect between downsteram region of inv and maf alignment')
    parser.add_argument('--intersect',metavar='<PDS>', help='File to save intersect between inv and maf alignment')
    parser.add_argument('--processMaf',metavar='<PDS>', help='File to save processed maf alignment')

    parser.set_defaults(entry_point=runFromArgs)

    args = parser.parse_args()
    sys.exit(args.entry_point(args))
