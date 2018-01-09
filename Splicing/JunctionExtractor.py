import argparse
from copy import deepcopy
import pysam
import sys

# Cigar operation dictory
CigarDict = {
    0: 'M', # BAM_CMATCH
    1: 'I', # BAM_CINS
    2: 'D', # BAM_CDEL
    3: 'N', # BAM_CREF_SKIP
    4: 'S', # BAM_CSOFT_CLIP
    5: 'H', # BAM_CHARD_CLIP
    6: 'P', # BAM_CPAD
    7: '=', # BAM_CEQUAL
    8: 'X', # BAM_CDIFF
    9: 'B'  # BAM_CBACK
}

class Junction:
    """ Class of a junction """
    def __init__(self, chrom, strand, start):
        """Init with chrom, chromStart and chromEnd"""
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = 0
        self.thick_start = self.start
        self.thick_end = self.end
        self.left_anchor = False
        self.right_anchor = False
        self.name = "Junc"
        self.read_count = 0
        self.score = "0"

    def __lt__(self, other):
        """ cmp function of junction """
        if self.chrom == other.chrom:
            if self.start == other.start:
                return self.end < other.end
            else:
                return self.start < other.start
        else:
            return self.chrom < other.chrom
    def __repr__(self):
        output = str(self.chrom)+'\t'+ \
                 str(self.thick_start)+'\t'+ \
                 str(self.thick_end)+'\t'+ \
                 str(self.name)+'\t'+ \
                 str(self.score)+'\t'+ \
                 str(self.strand)+'\t'+ \
                 str(self.thick_start)+'\t'+ \
                 str(self.thick_end)+'\t'+ \
                 "0,0,0\t2\t"+ \
                 str(self.start-self.thick_start)+','+ \
                 str(self.thick_end-self.end)+'\t'+ \
                 "0,"+str(self.end - self.thick_start)
        return output

    def filter(self, args):
        """ filter junctions """
        intron_len = self.end - self.start
        left_anchor = self.start - self.thick_start
        right_anchor = self.thick_end - self.end
        if intron_len < args.mini or intron_len > args.maxi:
            return False
        if left_anchor >= args.anchor:
            self.left_anchor = True
        if right_anchor >= args.anchor:
            self.right_anchor = True
        return True
def get_new_junction_name(SetLen):
    return 'Junction'+str(SetLen+1)

# Parser arguments
def parse_opt():
    parser = argparse.ArgumentParser(description='Extract Junction Reads from bamfile')
    parser.add_argument('-a', '--anchor', type=int, default=8, help='mininum anchor length')
    parser.add_argument('-b', '--bam', help='input bamfile', required=True)
    parser.add_argument('-i', '--mini', type=int, default=20, help='mininum intron length')
    parser.add_argument('-I', '--maxi', type=int, default=10000, help='maxinum intron length')
    parser.add_argument('-o', '--out', help='output file path')
    parser.add_argument('-r', '--region', help='junction reads from spcific regions')
    args = parser.parse_args()
    return args

# Get bamfile
def get_bamfile(bamfile):
    try:
        bam = pysam.AlignmentFile(bamfile, "rb")
        bam.check_index()
    except FileNotFoundError:
        msg='%s not found!' % bamfile
        print(msg)
    except AttributeError:
        msg='index of %s not found' % bamfile
    else:
        return bam

def get_junction_set():
    return [ v for v in sorted(JuncSet.values())]


# Add Junction
def add_junction(JuncSet, junction):
    junc = deepcopy(junction)
    #print "if same ",junc is junction
    #Check junction_qc
    if not junc.filter(args):
        return 
    tag = str(junc.chrom) + ':' + str(junc.start) + '-' + str(junc.end) + ':' + str(junc.strand)
    if tag in JuncSet:
        ori_junc = JuncSet[tag]
        junc.read_count = ori_junc.read_count + 1
        junc.score = str(junc.read_count)
        junc.name = ori_junc.name
        # Update thick_start
        if ori_junc.thick_start < junc.thick_start:
            junc.thick_start = ori_junc.thick_start
        if ori_junc.thick_end > junc.thick_end:
            junc.thick_end = ori_junc.thick_end
        # Update anchor information
        junc.left_anchor = junc.left_anchor or ori_junc.left_anchor
        junc.right_anchor = junc.right_anchor or ori_junc.right_anchor
        JuncSet[tag] = junc
    else:
        junc.name = get_new_junction_name(len(JuncSet))
        junc.read_count = 1
        junc.score = str(junc.read_count)
        JuncSet[tag] = junc
    #Add junction
# Parse alignment to junctions
def parse_junction_reads(JuncSet, chrom, ref_pos, strand, cigartuples):
    junc = Junction(chrom,strand,ref_pos)
    started_junction = False
    for cigar_oper, cigar_len in cigartuples:
        if CigarDict[cigar_oper] == 'M':
            if started_junction:
                junc.thick_end += cigar_len
            else:
                junc.start += cigar_len
        elif CigarDict[cigar_oper] == 'I':
            continue
        elif CigarDict[cigar_oper] == 'D':
            continue
        elif CigarDict[cigar_oper] == 'N':
            if started_junction:
                junc.thick_start = junc.end
                junc.start = junc.thick_end
                junc.end = junc.start + cigar_len
                junc.thick_end = junc.end
                #/For clarity - the next junction is now open
                started_junction = True
            else:
                junc.end = junc.start + cigar_len
                junc.thick_end = junc.end
                #Start the first one and remains started
                started_junction = True
        elif CigarDict[cigar_oper] == 'S':
            if started_junction:
                junc.start = junc.thick_end
                junc.thick_start = junc.start
            else:
                junc.thick_start = junc.start
            started_junction = False
        elif CigarDict[cigar_oper] == 'H':
            # NO hard clip is considered
            continue
        elif CigarDict[cigar_oper] == 'X':
            if started_junction:
                junc.start = junc.thick_end + cigar_len
                junc.thick_start = junc.start
            else:
                junc.start += cigar_len
                junc.thick_start = junc.start
            started_junction = False
        else:
            continue
    if started_junction:
        add_junction(JuncSet, junc)

# Get the strand based on 
def get_strand(aln):
    if aln.is_reverse:
        return '-'
    else:
        return '+'

# Get alignment 
def get_alignment(bam, JuncSet):
    alignment = bam.fetch()
    for aln in alignment:
        #print dir(aln)
        flag = aln.flag
        # flag filter 
        # read unmapped (0x4)
        # not primary alignment (0x100)
        # read is PCR or optical duplicate (0x400)
        if flag & 0x504:
            continue
        chrom = aln.reference_name
        ref_pos = aln.pos
        cigartuples = aln.cigartuples
        strand = get_strand(aln)
        # only one cigar in cigar string
        if len(cigartuples) == 1:
            continue
        
        parse_junction_reads(JuncSet, chrom, ref_pos, strand, cigartuples)

if __name__ == '__main__':
    args = parse_opt()
    bam = get_bamfile(args.bam)
    JuncSet = dict()
    get_alignment(bam, JuncSet)
    for junc in get_junction_set():
        print junc
