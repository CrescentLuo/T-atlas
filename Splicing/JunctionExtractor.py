import argparse
import pysam
import sys

# Parser arguments
parser = argparse.ArgumentParser(description='Extract Junction Reads from bamfile')
parser.add_argument('-a', '--anchor', help='mininum anchor length')
parser.add_argument('-i', '--mini', help='mininum intron length')
parser.add_argument('-I', '--maxi', help='maxinum intron length')
parser.add_argument('-o', '--out', help='output file path')
parser.add_argument('-r', '--region', help='junction reads from spcific regions')
args = parser.parse_args()

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
        
        


# Parse alignment to junctions
def parse_junction_reads(alignment):
    aln = alignment
    if aln.is_unmapped:
        return null
    reference = 
    ref_pos = 
    cigartuples = aln.cigartuples

    """123"""
    #int n_cigar = aln->core.n_cigar;
    if (n_cigar <= 1) // max one cigar operation exists(likely all matches)
        return 0;

    int chr_id = aln->core.tid;
    int read_pos = aln->core.pos;
    string chr(header->target_name[chr_id]);
    uint32_t *cigar = bam_get_cigar(aln);

    /*
    //Skip duplicates
    int flag = aln->core.flag;
    if(flag & 1024) {
        cerr << "Skipping read_pos " << read_pos << " flag " << flag << endl;
        return 0;
    }
    */

    Junction j1;
    j1.chrom = chr;
    j1.start = read_pos; //maintain start pos of junction
    j1.thick_start = read_pos;
    set_junction_strand(aln, j1);
    bool started_junction = false;

    started_junction = False
    for cigar_oper, cigar_len in cigar_tuple:
        if CigarDict[cigar_oper] == 'N':
            if !started_junction:
                
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    //Start the first one and remains started
                    started_junction = true;
                } else {
                    //Add the previous junction
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    j1.thick_start = j1.end;
                    j1.start = j1.thick_end;
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    //For clarity - the next junction is now open
                    started_junction = true;
                }
                break;
            if 
        elif CigarDict[cigar_oper] == 'M':
            if(!started_junction)
                    j1.start += len;
                else
                    j1.thick_end += len;
                break;
        elif CigarDict[cigar_oper] == 'X':
            if(!started_junction) {
                    j1.start += len;
                    j1.thick_start = j1.start;
                } else {
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    j1.start = j1.thick_end + len;
                    j1.thick_start = j1.start;
                }
                started_junction = false;
                break;
        elif CigarDict[cigar_oper] == 'S':
             if(!started_junction)
                    j1.thick_start = j1.start;
                else {
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    j1.start = j1.thick_end;
                    j1.thick_start = j1.start;
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
            add_junction(j1);
        } catch (const std::logic_error& e) {
            cout << e.what() << '\n';
        }
    }
    return 0;
}

# Get the strand based on 
def get_strand():
    return strand