import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='Length summary for minimum length.')
parser.add_argument('--summary', help='Length summary provided',
                    action="store_true")
parser.add_argument('--glob', help="Summary files or fastq files", type=str)
parser.add_argument('-m', '--max', help="Max length of the reads", type=int)
args = parser.parse_args()

file_set = glob.glob(args.glob)
if args.summary:
    for summary_file in file_set:
        with open(summary_file) as sum_file:
            read_sum = 0
            len_dict = {}
            prop_dict = {}
            for line in sum_file:
                read_cnt, read_len = line.rstrip().split()
                read_cnt = float(read_cnt)
                read_len = int(read_len)
                read_sum += read_cnt
                len_dict[read_len] = read_cnt
                prop_dict[read_len] = 0.0
            read_add = 0
            for read_len in range(args.max, -1, -1):
                read_add += len_dict[read_len]
                prop_dict[read_len] = read_add / read_sum
                if prop_dict[read_len] > 0.8:
                    print read_len
                    break
                
    


