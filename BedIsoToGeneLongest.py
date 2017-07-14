import sys
import argparse

parser = argparse.ArgumentParser(description='Trace gene to longest isoform bed format')
parser.add_argument('-g', '--gene', help='gene list file')
parser.add_argument('-b', '--bed', help='bed file of isoforms')
args = parser.parse_args()

with open(args.bed) as bed_file:
    for line in bed_file:
        sline = line.rstrip().split()
        