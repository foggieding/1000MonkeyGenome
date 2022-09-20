## Script: determine_inv.py
## Description: Determine the status of inversion.
## Author: Kevin Lee
## Date: 2022.01.12

# NOTE: Criteria used to determine inversion form was revised (2022-02-22)

import argparse
from utils import InputStream
from utils import PathType
import sys,os
import numpy as np
import re
import pickle
import operator

def warn(message):
    prog = os.path.basename(sys.argv[0])
    sys.stderr.write(prog + ": " + message + "\n")

def calculateSize(invName,invMaf):
    '''
    Calculate the size of forward and reverse alignment in chromosomes where inversion-intersected alignments span.
    :return: alnSize (dict)
    '''

    alnSize = dict()
    for aln in invMaf:
        chr2,strand,blocks = aln[2],aln[3],aln[4]
        chr2 = re.sub('^.*\.', '', chr2)
        sizeF=0;sizeR=0
        if strand == '+':
            sizeF = sum([b[2] for b in blocks])
        if strand == '-':
            sizeR = sum([b[2] for b in blocks])

        if chr2 in alnSize:
            alnSize[chr2] = tuple(map(operator.add,alnSize[chr2],(sizeF,sizeR)))
        else:
            alnSize[chr2] = (sizeF,sizeR)
    return alnSize

def determineChr(coord,index,alnSize):
    '''
    Determine alignment chromosome for inversion
     # Criteria:
     #  (1) The determined chromosome is the chromosome against which more than half the length of maf alignment in the interval mapped.
    :param coord:
    :param index:
    :param alnSize:
    :return: detChr
    '''
    invSize = invSize = coord[index][1] - coord[index][0]
    for chr, size in alnSize.items():
        if sum(size) / invSize >= 0.5:
            detChr = chr
            break
        else:
            continue
    try:
        detChr
    except NameError:
        detChr = None
    return detChr

def determineOri(detChr,alnSize,frac):
    '''
    Determine the orientation of inversion
    :param detChr:
    :param alnSize:
    :param frac:
    :return: strand
    '''
    sizeF,sizeR = alnSize[detChr]
    print(sizeF,sizeR,file=sys.stderr,flush=True)
    if sizeF / (sizeF+sizeR) >= frac:
        strand = "+"
    elif sizeR / (sizeF+sizeR) >= frac:
        strand = "-"
    else:
        strand = None
    return strand


def determineStatus(detChr, detChrUp, detChrDown, direction, directionUp, directionDown):
    '''
    Determine inversion status
    如果inversion跟上游或下游中任一区间染色体一致，并且方向为一正一反，则定义为inversion。
    :param detChr:
    :param detChrUp:
    :param detChrDown:
    :param direction:
    :param directionUp:
    :param directionDown:
    :return: status
    '''

    strandSet = ("+", "-")
    strandSet1 = ("+")
    strandSet2 = ("-")

    a = (detChr == detChrUp, detChr == detChrDown)
    b = (set([direction, directionUp]) == set(strandSet), set([direction, directionDown]) == set(strandSet))
    c = (set([direction, directionUp]) == set(strandSet1) or set([direction, directionUp]) == set(strandSet2),
         set([direction, directionDown]) == set(strandSet1) or set([direction, directionDown]) == set(strandSet2))

    if any(np.logical_and(a, b)):
        return "inverted"
    elif any(np.logical_and(a,c)):
        return "notInverted"
    else:
        return "uncertain"

def runFromArgs(args):
    print('Read inversion-related intervals',file=sys.stderr,flush=True)
    with open(args.inv,'rb') as f:
        inv = pickle.load(f)
    #print(inv)
    print('Read maf alignment intersected with inversion intervals',file=sys.stderr,flush=True)
    with open(args.maf,'rb') as f:
        maf = pickle.load(f)
    with open(args.mafUp,'rb') as f:
        mafUp = pickle.load(f)
    with open(args.mafDown,'rb') as f:
        mafDown = pickle.load(f)
    
    #print(maf)

    print("Determine inversion status",file=sys.stderr,flush=True)
    for invName,coord in inv:
        #if invName != "chr18_inv1": continue
        print("Processing {}".format(invName),file=sys.stderr,flush=True)
        invMaf = maf[invName]
        invMafUp = mafUp[invName]
        invMafDown = mafDown[invName]

        #print("inv size in hg38: ",invSize)
        #print("inv information: ",invName,coord,coord[index][1],coord[index][0])
        #print("maf alignment: ",maf[invName])

        # get the alignment length in forward and reverse direction in each chromosome that invertion-associated intervals span.
        alnSize = calculateSize(invName,invMaf)
        alnSizeUp = calculateSize(invName,invMafUp)
        alnSizeDown = calculateSize(invName,invMafDown)

        print("size in chromosomes: ", invName, alnSize,file=sys.stderr,flush=True)
        print("size of UP in chromosomes: ", invName, alnSizeUp,file=sys.stderr,flush=True)
        print("size of DOWN in chromosomes: ", invName, alnSizeDown,file=sys.stderr,flush=True)

        # determine chromosome
        detChr = determineChr(coord,0,alnSize)
        detChrUp = determineChr(coord,1,alnSizeUp)
        detChrDown = determineChr(coord,2,alnSizeDown)

        print("determinant chromosomes: ",detChr,file=sys.stderr,flush=True)
        print("determinant chromosomes UP: ",detChrUp,file=sys.stderr,flush=True)
        print("determinant chromosomes DOWN: ",detChrDown,file=sys.stderr,flush=True)

        # determine direction (orientation)
        direction = determineOri(detChr,alnSize,0.8) if detChr is not None else None
        directionUp = determineOri(detChrUp,alnSizeUp,0.8) if detChrUp is not None else None
        directionDown = determineOri(detChrDown,alnSizeDown,0.8) if detChrDown is not None else None

        print("determinant direction: ",direction,file=sys.stderr,flush=True)
        print("determinant direction UP: ",directionUp,file=sys.stderr,flush=True)
        print("determinant direction DOWN: ",directionDown,file=sys.stderr,flush=True)

        # determine inversion status
        status = determineStatus(detChr,detChrUp,detChrDown,direction,directionUp,directionDown)
        print("inversion status: ",status,file=sys.stderr,flush=True)

        # save to file
        l = [invName,detChr,detChrUp,detChrDown,direction,directionUp,directionDown,status]
        args.output.write('\t'.join([str(x) for x in l]) + '\n')

if __name__ == "__main__":
    descrition='Determine the status of inversion.'
    epilog='Please have your inversion-related intervals and maf alignment intersected with them ready.'
    parser = argparse.ArgumentParser(description=descrition,epilog=epilog)
    parser.add_argument('--inv',metavar='<PDS>',default=None,help='Inversion-related interval saved by pickle in get_invMaf_aln.py')
    parser.add_argument('--maf',metavar='<PDS>',default=None,help='MAF alignment intersected with inversion intervals, saved by pickle in get_invMaf_aln.py')
    parser.add_argument('--mafUp',metavar='<PDS>',default=None,help='MAF alignment intersected with region upstream inversion, saved by pickle in get_invMaf_aln.py')
    parser.add_argument('--mafDown',metavar='<PDS>',default=None,help='MAF alignment intersected with region downstream inversion, saved by pickle in get_invMaf_aln.py')
    parser.add_argument('--outDir',metavar='<DIR>',default='.', type=PathType(exists=True, type='dir'), help='Output directory')
    parser.add_argument('--output', metavar='<TXT>', type=argparse.FileType('w'), default=sys.stdout,
                       help='Output file to save inversion status (default: stdout)')

    parser.set_defaults(entry_point=runFromArgs)

    args = parser.parse_args()
    sys.exit(args.entry_point(args))
