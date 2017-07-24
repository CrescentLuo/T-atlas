import sys
import pysam
import argparse

parser = argparse.ArgumentParser(description="Splicing Index")
parser.add_argument('-i', '--input', help="input bam file")
parser.add_argument('-l', '--lnc', help="expressed lncRNAs")
parser.add_argument('-r', '--ref', help="bed file for referrence")
parser.add_argument('-d', '--dict', help="isoform gene dict")
args = parser.parse_args()

bam = pysam.AlignmentFile(args.input,"rb")

def calculateSI(iso,bam):
    chrom = iso[0]
    start = int(iso[1])
    end = int(iso[2])
    strand = iso[5]
    exonCnt = int(iso[9])
    exonLen = iso[10].rstrip(',').split(',')
    exonLen = [int(length) for length in exonLen]
    exonS = iso[11].rstrip(',').split(',')
    exonS = [int(s) for s in exonS]
    eeRead = 0.0
    eiRead = 0.0
    for i in range(1,exonCnt):
        exonPos = start + exonS[i]
        intronS = start + exonS[i-1] + exonLen[i-1]
        reads = bam.fetch(chrom, intronS - 100, exonPos + 3)
        for read in reads:
            if read.get_cigar_stats()[0][3] != 0:
                eeRead += 1
            else:
                eiRead += 1
    if eiRead != 0:
        splicing_index = eeRead / eiRead
    else:
        splicing_index = 0
    return splicing_index



iso_dict = dict()
with open(args.ref) as ref_file:
    for line in ref_file:
        sline = line.rstrip().split()
        iso_dict[sline[3]] = sline

gene_iso = dict()
with open(args.dict) as gene_iso_file:
    for line in gene_iso_file:
        sline = line.rstrip().split()
        exon_cnt = int(sline[8])
        exonS = sline[9].rstrip(',').split(',')
        exonS = [int(S) for S in exonS]
        exonE = sline[10].rstrip(',').split(',')
        exonE = [int(E) for E in exonS]
        iso_len = 0
        for i in range(exon_cnt):
            iso_len += (exonE[i] - exonS[i])
        if sline[0] not in iso_dict:
            gene_iso[sline[0]] = [sline[1],iso_len]
        else:
            if iso_len > iso_dict[sline[0]][1]:
                gene_iso[sline[0]] = [sline[1],iso_len] 
with open(args.lnc) as lncRNA_file:
    header = lncRNA_file.readline()
    print header.rstrip()
    for line in lncRNA_file:
        sline = line.rstrip().split()
        isoform = iso_dict[gene_iso[sline[0]][0]]
        si = calculateSI(isoform,bam)
        print line.rstrip()+'\t'+str(si)
        