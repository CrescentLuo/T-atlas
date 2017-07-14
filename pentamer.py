import sys
import argparse
import itertools
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
args = parser.parse_args()

def pentamerPermutation():
    nucbase = "ATCG"
    nucbase_pentamer_dict = dict()
    for idx,combination in enumerate(itertools.product(nucbase, repeat=5)):
        pentamer = ''.join(combination)
        nucbase_pentamer_dict[pentamer] = idx
    return nucbase_pentamer_dict

with open(args.bed) as bedFile:
    pentamer_dict = pentamerPermutation()
    fasta = BedTool(args.fasta)
    for line in bedFile:
        bedline = BedTool(line, from_string=True)
        get_fasta = bedline.sequence(fi=fasta, split=True, s=True)
        seq =  (open(get_fasta.seqfn).read()).split('\n')[1]
        #print seq
        seq = seq.upper()
        length = len(seq)
        pentamer_freq = [0.0] * (4 ** 5)
        isoform = line.rstrip().split()[4]
        for i in range(5, length):
            pentamer = seq[i-5:i]
            if pentamer in pentamer_freq:
                pentamer_freq[pentamer_dict[pentamer]] += 1
        for ind,cnt in enumerate(pentamer_freq):
            pentamer_freq[ind] = cnt / length * 1000 
        pentamer_freq = [str(freq) for freq in pentamer_freq]       
        print line.split()[3]+'\t'+'\t'.join(pentamer_freq)
        