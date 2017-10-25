import sys
import argparse

parser = argparse.ArgumentParser(description="Join two files together")
parser.add_argument('-a', help="file one")
parser.add_argument('-b', help="file two")
parser.add_argument('-j', type=int, help="join field")
args = parser.parse_args()

with open(args.a) as left_file, open(args.b) as right_file:
    left_dict = dict()
    for line in left_file:
        sline = line.rstrip().split()
        left_dict[sline[args.j - 1]] = line.rstrip()
    for line in right_file:
        sline = line.rstrip().split()
        if sline[args.j - 1] in left_dict:
            print left_dict[sline[args.j - 1]]+'\t'+line.rstrip()
