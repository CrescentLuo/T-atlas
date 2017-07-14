import sys
import argparse
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
args = parser.parse_args()


with open(args.bed) as bedFile:
    fasta = BedTool(args.fasta)
    for line in bedFile:
        bedline = BedTool(line, from_string=True)
        get_fasta = bedline.sequence(fi=fasta, split=True, s=True)
        seq =  (open(get_fasta.seqfn).read()).split('\n')[1]
        #print seq
        seq = seq.upper()
        countA = float(seq.count('A'))
        countT = float(seq.count('T'))
        countC = float(seq.count('C'))
        countG = float(seq.count('G'))
        seq_len= len(seq)
        percentageA = countA / seq_len 
        percentageT = countT / seq_len
        percentageC = countC / seq_len
        percentageG = countG / seq_len
        print line.split()[3],'\t',percentageA,'\t',percentageT,'\t',percentageC,'\t',percentageG
        