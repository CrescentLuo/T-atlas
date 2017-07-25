import sys
import argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='get fasta from bed files and count ATCG percentage')
parser.add_argument('-b', '--bed', required=True, help='input bed file')
parser.add_argument('-f', '--fasta', required=True, help='fasta file')
parser.add_argument('-k', '--repeat', required=True, help='set K for kmer')
parser.add_argument('-s', '--split', default=True, help='split to get mature sequence')
parser.add_argument('-p', '--process', default=5, help='number of processes')
parser.add_argument('--overlap', action='store_true', help='Count kmer with or without overlapping')
args = parser.parse_args()

def kmerPermutation(K):
    nucbase = "ATCG"
    nucbase_kmer_dict = dict()
    for idx,combination in enumerate(itertools.product(nucbase, repeat=K)):
        kmer = ''.join(combination)
        nucbase_kmer_dict[kmer] = idx
    header = "isoform\t"+'\t'.join([ Kmer for Kmer in nucbase_kmer_dict])
    print header
    return nucbase_kmer_dict

def kmerFreq( K, line):
    sline = line.rstrip().split()
    chrom = sline[0]
    start = sline[1]
    end = sline[2]
    exonCnt = int(sline[10])
    exonlen = sline[10].rstrip(',').split(',')
    exonlen = [int(length) for length in exonlen]
    exonS = sline[11].rstrip(',').split(',')
    exonS = [int(s) for s in exonS]
    strand = sline[5]
    gene_seq = records[chrom].seq[int(sline[1]):int(sline[2])].upper()
<<<<<<< HEAD
=======
    gene_seq_str = str(gene_seq)
>>>>>>> 3037e51f0e1bdad400d9c06b320d97051dc8cada
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
<<<<<<< HEAD
    isoform = line.rstrip().split()[4]
    for kmer in kmer_dict:
        kmer_freq[kmer_dict[kmer]] += spliced_seq
=======
    isoform = sline[4]
    for kmer in kmer_dict:
        if args.overlap:
            kmer_ferq[kmer_dict[kmer]] = spliced_seq.count_overlap(kmer) + 0.0
        else:
            kmer_freq[kmer_dict[kmer]] = spliced_seq.count(kmer) + 0.0
>>>>>>> 3037e51f0e1bdad400d9c06b320d97051dc8cada
    for ind,cnt in enumerate(kmer_freq):
        kmer_freq[ind] = cnt / spliced_length * 1000 
    kmer_freq = [str(freq) for freq in kmer_freq]  
    print line.split()[3]+'\t'+'\t'.join(kmer_freq)
    #return line.split()[3]+'\t'+'\t'.join(kmer_freq)
records = SeqIO.to_dict(SeqIO.parse(open(args.fasta), 'fasta'))
with open(args.bed) as bedFile:
    K = int(args.repeat)
    kmer_dict = kmerPermutation(K)
    pool = Pool(processes = int(args.process))
    for line in bedFile:
        pool.apply_async(kmerFreq, args=(K, line,))
    pool.close()
    pool.join()
