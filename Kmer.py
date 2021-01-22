__author__="Zheng Luo"
import sys
import argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
parser.add_argument('-k', '--repeat', required=True, type=int, help='set K for kmer')
parser.add_argument('-s', '--split', default=True, help='split to get mature sequence')
parser.add_argument('-p', '--process', default=5, type=int,help='number of processes')
parser.add_argument('-c', '--count', default=False, action='store_true', help='use count instead of frequency')
parser.add_argument('-o', '--output', default="kmer_freq.txt", help='output file path')
parser.add_argument('--overlap', action='store_true', help='Count kmer with or without overlapping')
args = parser.parse_args()

def kmerPermutation(K):
    nucbase = "ATCG"
    nucbase_kmer_dict = dict()
    header = "isoform"
    for idx,combination in enumerate(itertools.product(nucbase, repeat=K)):
        kmer = ''.join(combination)
        nucbase_kmer_dict[kmer] = idx
        header = header + '\t' + kmer
    return header+'\n', nucbase_kmer_dict

def kmerFreq(isoform):
    K = args.repeat
    sline = isoform.rstrip().split()
    chrom = sline[0]
    start = sline[1]
    end = sline[2]
    exonCnt = int(sline[9])
    exonlen = sline[10].rstrip(',').split(',')
    exonlen = [int(length) for length in exonlen]
    exonS = sline[11].rstrip(',').split(',')
    exonS = [int(s) for s in exonS]
    strand = sline[5]
    gene_seq = records[chrom].seq[int(sline[1]):int(sline[2])].upper()
    gene_seq_str = str(gene_seq)
    spliced_seq = ""
    for i in range(exonCnt):
        spliced_seq = spliced_seq + gene_seq_str[exonS[i]:(exonS[i] + exonlen[i])]
    if strand == '-':
        gene_seq = gene_seq.reverse_complement()
        spliced_seq = Seq(spliced_seq).reverse_complement()
    else:
        spliced_seq = Seq(spliced_seq)
    spliced_length = len(spliced_seq)
    gene_length = len(gene_seq)
    kmer_freq = [0.0] * (4 ** K)
    iso = sline[3]
    for kmer in kmer_dict:
        if args.overlap:
            kmer_freq[kmer_dict[kmer]] = spliced_seq.count_overlap(kmer) + 0.0
        else:
            kmer_freq[kmer_dict[kmer]] = spliced_seq.count(kmer) + 0.0
    if not(args.count):
        for ind,cnt in enumerate(kmer_freq):
            kmer_freq[ind] = cnt / spliced_length * 1000 
    kmer_freq = [str(freq) for freq in kmer_freq]  
    return iso+'\t'+'\t'.join(kmer_freq)+'\n'
    #return line.split()[3]+'\t'+'\t'.join(kmer_freq)
records = SeqIO.to_dict(SeqIO.parse(open(args.fasta), 'fasta'))

with open(args.bed) as bedFile, open(args.output,"w") as output:
    isoforms = bedFile.readlines()
    K = int(args.repeat)
    header, kmer_dict = kmerPermutation(K)
    output.write(header)
    pool = Pool(processes = int(args.process))
    for freq in pool.imap(kmerFreq, isoforms):
        output.write(freq)
    pool.close()
    pool.join()
