import argparse
import numpy as np
import sys


parser = argparse.ArgumentParser(description='rMATs result filter')
parser.add_argument('-i', '--input', required=True, help='rMATs results as input file')
parser.add_argument('-c', '--count', type=float, default=0, help='mininum read count for AS events')
parser.add_argument('--pval', type=float, default=0.05, help='pval cutoff for significant AS events')
parser.add_argument('--FDR', type=float, default=0.5, help='FDR cutoff')
parser.add_argument('--psi', type=float, default=0.0, help='mininum delta psi')
args = parser.parse_args()

with open(args.input) as as_events:
    for line in as_events:
        event = line.rstrip().split()
        ijc_samp1 = [float(c) for c in event[12].split(',')]
        sjc_samp1 = [float(c) for c in event[13].split(',')]
        ijc_samp2 = [float(c) for c in event[14].split(',')]
        sjc_samp2 = [float(c) for c in event[15].split(',')] 
        ijc_samp1_flag = False
        sjc_samp1_flag = False
        ijc_samp2_flag = False
        sjc_samp2_flag = False

        if np.mean(ijc_samp1) >= args.count:
                ijc_samp1_flag = True
        if np.mean(ijc_samp2) >= args.count:
                ijc_samp2_flag = True
        if np.mean(sjc_samp1) >= args.count:
                sjc_samp1_flag = True
        if np.mean(sjc_samp2) >= args.count:
                sjc_samp2_flag = True
        if not(ijc_samp1_flag or sjc_samp1_flag) or not(ijc_samp2_flag or sjc_samp2_flag):
            continue

        pval = float(event[18])
        if pval > args.pval:
            continue
        fdr = float(event[19])
        if fdr > args.FDR:
            continue
        psi = abs(float(event[22]))

        if psi < args.psi:
            continue
        print line.rstrip()