import sys
import argparse

parser = argparse.ArgumentParser(description="Extract exons and introns from bed file")
parser.add_argument('-i', '--input', help="input bed file")
parser.add_argument('-n', '--name', help="out put name prefix")
args = parser.parse_args()

with open(args.input) as input_bed, \
     open((args.name + "_exon.bed"),"w") as exon_file, \
     open((args.name + "_intron.bed"),"w") as intron_file:  

    for line in input_bed:
        chrom, start, end, name, score, strand, cdssta, cdsend, rgb, exonCnt, sizes, offsets = line.rstrip().split()
        start = int(start)
        end = int(end)
        exonCnt = int(exonCnt)
        sizes = [int(size) for size in sizes.rstrip(',').split(',')]
        offsets = [int(offset) for offset in offsets.rstrip(',').split(',')]
        for i in xrange(int(exonCnt)):
            regionS = start + offsets[i]
            regionE = start + offsets[i] + sizes[i]
            exon_name = name+'_exon_'+str(i+1)
            exon_size = str(sizes[i])+','
            exon = [chrom, regionS, regionE, exon_name, score, strand, regionS, regionE, rgb, "1", exon_size, '0,']
            exon = ('\t'.join([str(si) for si in exon]) + '\n') 
            exon_file.write(exon)
        for i in xrange(int(exonCnt - 1)):
            regionS = start + offsets[i]+sizes[i]
            regionE = start + offsets[i+1]
            intron_name = name+'_intron_'+str(i+1)
            intron_size = offsets[i+1] - sizes[i] - offsets[i]
            if intron_size == 0:
                continue
            intron_size = str(offsets[i+1] - sizes[i] - offsets[i])+','
            intron =[chrom, regionS, regionE, intron_name, score, strand, regionS, regionE, rgb, "1", intron_size,'0,']
            intron = ('\t'.join([str(si) for si in intron]) + '\n') 
            intron_file.write(intron)