import sys
import argparse
import itertools

parser = argparse.ArgumentParser(description="melt pentamer file to top pentamers")
parser.add_argument('-i', '--input', help='input pentamer file')

args = parser.parse_args()
def pentamerPermutation():
    nucbase = "ATCG"
    nucbase_pentamer_dict = dict()
    for idx,combination in enumerate(itertools.product(nucbase, repeat=5)):
        pentamer = ''.join(combination)
        nucbase_pentamer_dict[idx] = pentamer
    return nucbase_pentamer_dict

with open(args.input) as input_file:
    pentamer_loc = pentamerPermutation()
    for line in input_file:
        sline = line.rstrip().split()
        gene_id = sline[0]
        symbol = sline[1]
        wt_nuc_prop = sline[25]
        fast_nuc_prop = sline[15]
        pentamer_rank = list()
        for i in range(39,1063):
            pentamer_rank.append([float(sline[i]),i-39,0])
        pentamer_rank = sorted(pentamer_rank,reverse=True)
        for i in range(1024):
            print gene_id+'\t'+symbol+'\t'+wt_nuc_prop+'\t'+fast_nuc_prop+'\t'+pentamer_loc[pentamer_rank[i][1]]+'\t'+str(i+1)+'\t'+str(pentamer_rank[i][0])