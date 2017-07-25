import sys
import argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
parser.add_argument('-k', '--repeat', required=True, help='set K for kmer')
parser.add_argument('-s', '--split', default=True, help='split to get mature sequence')
parser.add_argument('-p', '--process', default=5, help='number of processes')
parser.add_argument('-o', '--overlaping', default=True, help='Count kmer with or without overlapping')
args = parser.parse_args()

def kmerPermutation(K):
    nucbase = "ATCG"
    nucbase_kmer_dict = dict()
    for idx,combination in enumerate(itertools.product(nucbase, repeat=K)):
        kmer = ''.join(combination)
        nucbase_kmer_dict[kmer] = idx
    return nucbase_kmer_dict

def kmerFreq( K, line):
    sline = line.rstrip().split()
    chrom = sline[0]
    start = sline[1]
    end = sline[2]
    exonCnt = int(sline[9])
    exonlen = sline[10].rstrip(',').split(',')
    exonS = sline[11].rstrip(',').split(',')
    strand = sline[5]
    gene_seq = records[chrom].seq[int(sline[1]):int(sline[2])].upper()
    spliced_seq = ""
    for i in range(exonCnt):
        spliced_seq = spliced_seq + gene_seq[exonS[i]:(exonS[i] + exonlen[i])]
    print len(spliced_seq),'\t',spliced_seq
    seq =  (open(get_fasta.seqfn).read()).split('\n')[1]
    #print seq
    seq = seq.upper()
    #print seq
    length = len(seq)
    kmer_freq = [0.0] * (4 ** K)
    isoform = line.rstrip().split()[4]
    for kmer in kmer_dict:
        kmer_freq[kmer_dict[kmer]] += spliced_seq
    for ind,cnt in enumerate(kmer_freq):
        kmer_freq[ind] = cnt / length * 1000 
    kmer_freq = [str(freq) for freq in kmer_freq]  
    print line.split()[3]+'\t'+'\t'.join(kmer_freq)
    #return line.split()[3]+'\t'+'\t'.join(kmer_freq)
records = SeqIO.to_dict(SeqIO.parse(open(args.fasta), 'fasta'))
with open(args.bed) as bedFile:
    K = int(args.repeat)
    kmer_dict = kmerPermutation(K)
    for line in bedFile:
        kmerFreq(K, line)
        