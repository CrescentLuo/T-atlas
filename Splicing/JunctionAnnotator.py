import argparse
from copy import deepcopy
import pysam
import sys

def parse_opt():
    parser = argparse.ArgumentParser(description='Extract Junction Reads from bamfile')
    parser.add_argument('-a', '--anchor', type=int, default=8, help='mininum anchor length')
    parser.add_argument('-b', '--bam', help='input bamfile', required=True)
    parser.add_argument('-i', '--mini', type=int, default=70, help='mininum intron length')
    parser.add_argument('-I', '--maxi', type=int, default=500000, help='maxinum intron length')
    parser.add_argument('-o', '--out', help='output file path')
    parser.add_argument('-r', '--region', help='junction reads from spcific regions')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_opt()