import argparse
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

class Junction():
    """ Class of a junction """
    def __init__(self, chrom, start, end):
        """Init with chrom, chromStart and chromEnd"""
        self.chrom = chrom
        self.start = start
        self.end = end

JuncSet = dict()

# Parser arguments
def parse_opt():
    parser = argparse.ArgumentParser(description='Extract Junction Reads from bamfile')
    parser.add_argument('-a', '--anchor', help='mininum anchor length')
    parser.add_argument('-i', '--mini', help='mininum intron length')
    parser.add_argument('-I', '--maxi', help='maxinum intron length')
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

# Add Junction
def add_junction(junc):
    #Check junction_qc
    """
    if(!junction_qc(junc)) {
        return 0;
    }
    """
    tag = junc.chrom + ':' + str(junc.start) + '-' + str(junc.end) + ':' + junc.strand
    if tag not in JunctSet:
        junc.name = get_new_junction_name(JunctSet)
        junc.read_count = 1
        junc.score = str(junc.read_count)
    else:
        ori_junc = JuncSet[tag]
        junc.read_count = ori_junc.read_count + 1
        junc.score = str(junc.read_count)
        junc.name = ori_junc.name
        # Update thick_start
        if ori_junc.thick_start < junc.thick_start:
            junc.thick_start = ori_junc.thick_start;
        if ori_junc.thick_end > junc.thick_end:
            junc.thick_end = ori_junc.thick_end
        # Update anchor information
        junc.left_anchor = junc.left_anchor or ori_junc.left_anchor
        junc.right_anchor = junc.right_anchor or ori_junc.right_anchor
    #Add junction 
    JuncSet[tag] = junc

# Parse alignment to junctions
def parse_junction_reads(chrom, ref_pos, strand, cigartuples):
    junc = Junction(chrom,strand,ref_pos,ref_pos)
    started_junction = False
    for cigar_oper, cigar_len in cigar_tuple:
        if CigarDict[cigar_oper] == 'N':
            if not started_junction:
                junc.end = junc.start + cigar_len
                junc.thick_end = junc.end
                #Start the first one and remains started
                started_junction = true;
            else:
                #Add the previous junction
                try：
                    add_junction(junc)
                except (const std::logic_error& e
                        cout << e.what() << '\n'
                junc.thick_start = junc.end;
                junc.start = junc.thick_end;
                junc.end = junc.start + cigar_len;
                junc.thick_end = junc.end;
                #/For clarity - the next junction is now open
                started_junction = true
        elif CigarDict[cigar_oper] == 'M':
            if not started_junction:
                junc.start += cigar_len
            else:
                junc.thick_end += cigar_len
        elif CigarDict[cigar_oper] == 'X':
            if not started_junction:
                junc.start += cigar_len
                junc.thick_start = junc.start
            else:
                    try {
                        add_junction(junc);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    junc.start = junc.thick_end + len;
                    junc.thick_start = junc.start;
                }
                started_junction = false;
                break;
        elif CigarDict[cigar_oper] == 'S':
             if not started_junction
                    junc.thick_start = junc.start;
                else {
                    try {
                        add_junction(junc);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    junc.start = junc.thick_end;
                    junc.thick_start = junc.start;
                }
                started_junction = false;
                break;
        elif CigarDict[cigar_oper] == 'H':
            # NO hard clip is considered
            break
        else:
            "Unknown Cigar"
        if started_junction:
            try {
            add_junction(junc);
        } catch (const std::logic_error& e) {
            cout << e.what() << '\n';
        }
    }
    return 0;
}

# Get alignment 
def get_alignment(bam):
    alignment = bam.fetch()
    for aln in alignment:
        flag = aln.flag
        
        # flag filter 
        # read unmapped (0x4)
        # not primary alignment (0x100)
        # read is PCR or optical duplicate (0x400)
        
        if flag & 0x504:
            continue

        chrom = aln.get_reference_name()
        ref_pos = aln.get_pos()
        cigartuples = aln.cigartuples
        strand = get_strand(aln)
        # only one cigar in cigar string
        if len(cigartuples) == 1:
            continue
        
        parse_junction_reads(chrom, ref_pos, strand, cigartuples)
        
# Get the strand based on 
def get_strand(aln):
    if aln.is_reverse:
        return '-'
    else:
        return '+'

if __name__ == '__main__':
    args = parse_opt()
    bam = get_bamfile(args.bamfile)
    

