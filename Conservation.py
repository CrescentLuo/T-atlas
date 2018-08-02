import sys
import argparse
import pyBigWig as pyBW
import numpy as np

def arg_parser():
    parser = argparse.ArgumentParser(description="Calculation of PhastCons of genes using bigwig file")
    parser.add_argument('-b', '--bw', required=True, help='Input bigwig file')
    parser.add_argument('-a', '--anno', required=True, help='Annoation file(Bed format) of gene genomic location')
    args = parser.parse_args()
    return args


def anno_parser(anno_file):
    anno_dict = dict()
    anno_maxLen = dict()
    with open(anno_file) as input_file:
        for line in input_file:
            sline = line.rstrip().split()
            chrom = sline[0]
            chromS = int(sline[1])
            chromE = int(sline[2])
            geneID = sline[3]
            exonCnt = int(sline[9])
            exonSizes = [int(s) for s in sline[10].rstrip(',').split(',')]
            eoxnStarts = [int(s) for s in sline[11].rstrip(',').split(',')]
            exonLen = 0
            exonSets = list()
            for i in range(exonCnt):
                exonLen += exonSizes[i]
                exonSets.append((chromS+eoxnStarts[i], chromS+eoxnStarts[i]+exonSizes[i]))
            if geneID in anno_dict:
                if exonLen > anno_maxLen[geneID]:
                    anno_dict[geneID] = [chrom, exonSets]
                    anno_maxLen[geneID] = exonLen
            else:
                anno_dict[geneID] = [chrom, exonSets]
                anno_maxLen[geneID] = exonLen
    return anno_dict


if __name__ == '__main__':
    args = arg_parser()
    anno_dict = anno_parser(args.anno)
    with pyBW.open(args.bw) as bw:
        for gene in anno_dict:
            chrom = anno_dict[gene][0]
            exonSets = anno_dict[gene][1]
            exon_phastcons_sets = list()
            if chrom not in bw.chroms():
                continue
            for s,e in exonSets:
                phastcons_score = bw.stats(chrom, s, e)[0]
                if phastcons_score:
                    exon_phastcons_sets.append(phastcons_score)
            print gene, np.mean(exon_phastcons_sets)
        

