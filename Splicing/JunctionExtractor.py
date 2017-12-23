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
    
# Get the strand based on 
def get_strand():
    return strand