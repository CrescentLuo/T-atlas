import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='trans relative position to chromosomal coordinate')
    parser.add_argument('-b', '--bed', help='GENCODE annotation in bed format')
    parser.add_argument('-i', '--input', help='Input region file')
    return parser.parse_args()

def load_annotation(anno_file):
    """
    GENCODE annotation in bed format:
    "chr1    360056  366052  ENST00000635159.1       0       +       366052  366052  0       2       112,882,        0,5114,"
    """
    anno_set = dict()
    with open(anno_file) as anno:
        for line in anno:
            sline = line.rstrip().split()
            anno_set[sline[3]] = sline
    return anno_set

def coor_trans(input_file, anno_set):
    """
    Input file format (first 9 columns):
    "ENST00000508424.5  482-581 0.8807275126    \
    ACATAGCTTAGGTAGTGATTTTTGTTCATAAGAGTATCTCATAGTTTAATTATTCTGTCTATTTAAGAAATTACCTGGGAGTTAAAATTGGTCTTGGTAT \
    SEMA6A-AS1  chr5    antisense   ENST00000508424.5   ENSG00000248445.5 ..."
    """
    with open(input_file) as input_data:
        for line in input_data:
            sline = line.rstrip().split()
            if sline[0] not in anno_set:
                continue
            else:
                anno = anno_set[sline[0]]
                chrom = anno[0]
                chrom_s = int(anno[1])
                chrom_e = int(anno[2])
                
                
                

if __name__ == '__main__':
    args = parse_args()
    anno_set = load_annotation(args.bed)
    
