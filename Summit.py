#!/usr/bin/env python
import sys 
import pysam
import argparse
from collections import defaultdict
from multiprocessing import Pool

def argumentparse():
    parser = argparse.ArgumentParser(description="Peak finder")
    parser.add_argument('-c', '--control', help="control bam file")
    parser.add_argument('-t', '--target', help="target bam file")
    parser.add_argument('-s', '--size', help="sliding window size")
    parser.add_argument('-g', '--gtf', help="gene annotation file")
    parser.add_argument('-b', '--bed', help="bed format annotation file")
    parser.add_argument('-p', '--process', type=int, help="number of process ")
    parser.add_argument('-o', '--outfile', default="output.txt", help="output file")
    parser.add_argument('-e', '--exp', type=float, default=5, help="expression cutoff")
    args = parser.parse_args()
    return args

def parsebed(line):
    items = line.split()
    if len(items) == 12:
        chrom, start, end, name, score, strand, cdssta, cdsend, rgb, exonCnt, sizes, offsets = item
        start = int(start)
        end = int(end)
        exonCnt = int(exonCnt)
        sizes = [int(size) for size in sizes.split(',')]
        offsets = [int(offset) for offset in offsets.split(',')]
        regions = [chrom]
        for i in xrange(exonCnt):
            regions.append([start + offsets[i], start + offsets[i] + sizes[i]])
        introns = []
        for i in xrange(exonCnt - 1):
            introns.append([chrom, [start + offsets[i] + sizes[0], start + offsets[i + 1]], i + 1])
    else:
        chrom, start, end, name, score, strand = items
        regions = [chrom, [int(start), int(end)]]
        introns = []
    return (name, strand, regions, chrom, introns)

def cal_p(regions, exp_lvl, strand, name, size, repeats, FDR, t_bam, c_bam):
    for region in regions:
        mapped_reads = t_bam.fetch(chrom, region[0], region[1])
        region_q = list()
        step_q = list()
        base_win_set = list()
        win_start = region[0]
        win_end = region[0] + winSize
        win_cnt = 0.0
        extending = False
        for read in mapped_reads:
            if win_end >= region[1]:
                break
            if read.pos >= win_start and read.pos < win_end:
                if extending:
                    step_q.append(read.pos)
                else:
                    q.put(read.pos)
            else:
                win_extend_cnt = q.qsize()
                if extending:
                    win_exp = win_cnt / ((win_end - win_start) * t_bam.mapped_reads) * (10**9)
                else:
                    win_exp = win_cnt / ((win_end - win_start - step) * t_bam.mapped_reads) * (10**9)
                # Count the expression level of window
                # If window expression is higher than cutoff then open extending frame
                if win_exp >= exp_lvl:
                    # If the window is extended then calculate current extending frame 
                    # Decide whether to extend one more step
                    if extending:
                        step_exp = step_cnt / (step * t_bam.mapped_reads) * (10**9)
                        if win_exp * 1.5 >= step_exp and win_exp * 0.667 <= step_exp:
                            win_end = win_end + step_exp
                            extending =  True
                        else:
                            base_win_set.append(chrom, win_start, win_end, win_exp)
                            extending = False
                    else:
                        win_end = win_end + step
                        extending = True
                else:
                    extending = False
                    win_start += step
                    win_end += step
                    while q[0].pos < win_start:
                        q.pop()
        if extending:
            base_win_set.append(chrom, win_start, win_end, win_exp)
        p = list()
        win_set = list()
        for i in xrange(len(base_win_set)):
            if i == 0:
                p = base_win_set[i]
            else:
                if base_win_set[i][1] - p[i-1][2]:
                    p[i][2] = base_win_set[i][1]
                    if i == len(base_win_set - 1):
                        win_set.append(p)
                else:
                    win_set.append(p)
                    if i == len(base_win_set - 1):
                        win_set.append(base_win_set[i])
    return win_set
                
def main():
    t_bam = pysam.AlignmentFile(args.target,"rb")
    c_bam = pysam.AlignmentFile(args.control,"rb")
    outfile = args.outfile
    exp_cutoff = args.exp
    with open(args.bed) as ref_file:
        p = Pool(processes=args.process)
        results = []
        for line in ref_file:
            name, strand, regions, chrom, introns = parsebed(line.rstrip())
            result.append(p.apply_async(cal_p, args=(regions, exp_cutoff, strand, name, size, repeats, FDR, t_bam, c_bam)))
            for intron in introns:
                intron_name = '%s_intron%i' % (name, intron[-1])
                result.append(p.apply_async(cal_p, args=(intron[:-1], exp_cutoff, strand, intron_name, size, repeats, FDR, t_bam, c_bam)))
        p.close()
        p.join()
        for res in results:
            for win in res.get()
                print win
            