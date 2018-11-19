import sys
import argparse

def arg_parser():
    parser = argparse.ArgumentParser(description='Trace gene to longest isoform bed format')
    parser.add_argument('-g', '--gene', help='gene list file')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gp', help='Use genePred annotation as input file')
    group.add_argument('--bed', help='Use bed annotation as input file')
    args = parser.parse_args()
    return args

def open_bed(bedfile):
    with open(bedfile) as bed:
        isoform_dict = dict()
        for line in bed:
            sline = line.rstrip().split()
            gene = sline[3]
            exon_cnt = int(sline[9])
            exon_sizes = [int(s) for s in sline[10].rstrip(',').split(',')]
            exons_starts = [int(e) for e in sline[11].rstrip(',').split(',')]
            exon_len = 0
            for i in range(exon_cnt):
                exon_len += exon_sizes[i]
            if gene in isoform_dict:
                if exon_len >= isoform_dict[gene][1]:
                    isoform_dict[gene] = (line.rstrip(),exon_len)
            else:
                isoform_dict[gene] = (line.rstrip(),exon_len)
        return isoform_dict

def open_gp(genePredfile):
    with open(genePredfile) as gp:
        isoform_dict = dict()
        for line in gp:
            sline = line.rstrip().split()
            iso = sline[0]
            gene = sline[11]
            exon_cnt = int(sline[7])
            exons_starts = [int(s) for s in sline[8].rstrip(',').split(',')]
            exons_ends = [int(e) for e in sline[9].rstrip(',').split(',')]
            exon_len = 0
            for i in range(exon_cnt):
                exon_len += exons_ends[i] - exons_starts[i]
            if iso in isoform_dict:
                if exon_len >= isoform_dict[gene][1]:
                    isoform_dict[gene] = (iso,exon_len)
            else:
                isoform_dict[gene] = (iso,exon_len)
    return isoform_dict
                


if __name__ == '__main__':
    args = arg_parser()
    if args.gp:
        isoform_dict = open_gp(args.gp)
        for gene in isoform_dict:
            print gene,'\t',isoform_dict[gene][0]
    else:
        isoform_dict = open_bed(args.bed)
        for gene in isoform_dict:
            print isoform_dict[gene][0]

    

