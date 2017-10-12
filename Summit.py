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
            t_bam.fetch()
        
        

def main():
    t_bam = pysam.AlignmentFile(args.target,"rb")
    c_bam = pysam.AlignmentFile(args.control,"rb")
    outfile = args.outfil
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
            if res.get()
                outf.write(r.get())
            