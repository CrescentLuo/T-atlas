import sys
import argparse

parser = argparse.ArgumentParser(description="Merging iCLIP binding sites into clusters")
parser.add_argument('-s', '--site', help="position of binding sites")
parser.add_argument('-r', '--range' help="Genomic region of bingding region:'chr1:1000-2000'")
parser.add_argument('--min', type=int, default=50, help="minimus distance between two sites")
args = parser.parse_args()

bs_set = list[]
with open(args.site) as site_file：
    r_chrom = args.r.split(':')[0]
    r_start = int(args.r.split(':')[1].split('-')[0])
    r_end = int(args.r.split(':')[1].split('-')[1])
    for line in site_file:
        chrom, start, end, ts, cnt, strand, pval = line.rstrip().split()
        start = int(start)
        end = int(end)
        cnt = int(cnt)
        pval = float(pval)
        if chrom != r_chrom or start < r_start or start > r_end:
            continue
        else:
            bs_set.append([chrom, start, end, ts, cnt, strand, pval])
    for i in xrange(len(bs_set)):
        if i == 0:
            p = bs_set[i]
        if i > 0:
            q = bs_set[i]
            if q[1] - p[2] <= args.min:
                p[2] = q[2]
                p[4] += q[4]
            else:
                p = [str(info) for info in p]
                print '\t'.join(p[0:3]) + '\t' + '\t'.join(p[4:6])
                p = q
                if i ==  len(bs_set) - 1:
                    print '\t'.join(q[0:3]) + '\t' + '\t'.join(q[4:6]) 
        
