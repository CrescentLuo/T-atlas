import sys
import argparse
import itertools
import pybedtools
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
parser.add_argument('-k', '--repeat', required=True, help='set K for kmer')
parser.add_argument('-s', '--split', default=True, help='split to get mature sequence')
parser.add_argument('-p', '--process', default=5, help='number of processes')
args = parser.parse_args()

def kmerPermutation(K):
    nucbase = "ATCG"
    nucbase_kmer_dict = dict()
    for idx,combination in enumerate(itertools.product(nucbase, repeat=K)):
        kmer = ''.join(combination)
        nucbase_kmer_dict[kmer] = idx
    return nucbase_kmer_dict

def kmerFreq( K, line):
    bedline = BedTool(line, from_string=True)
    fasta = pybedtools.example_filename(args.fasta)
    get_fasta = bedline.sequence(fi=fasta, split=True, s=True)
    #print get_fasta
    #print get_fasta
    seq =  (open(get_fasta.seqfn).read()).split('\n')[1]
    #print seq
    seq = seq.upper()
    #print seq
    length = len(seq)
    kmer_freq = [0.0] * (4 ** K)
    isoform = line.rstrip().split()[4]
    for i in range(K, length):
        kmer = seq[i-K:i]
        if kmer in kmer_dict:
            kmer_freq[kmer_dict[kmer]] += 1
    for ind,cnt in enumerate(kmer_freq):
        kmer_freq[ind] = cnt / length * 1000 
    kmer_freq = [str(freq) for freq in kmer_freq]  
    print line.split()[3]+'\t'+'\t'.join(kmer_freq)
    #return line.split()[3]+'\t'+'\t'.join(kmer_freq)

with open(args.bed) as bedFile:
    pybedtools.set_tempdir('./')
    K = int(args.repeat)
    kmer_dict = kmerPermutation(K)
    for line in bedFile:
        kmerFreq(K, line)
        