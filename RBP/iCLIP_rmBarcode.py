import sys
import argparse
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description='Remove and record barcode')
    parser.add_argument('-p', '--pattern', required=True,
                        help='pattern of barcode"XXXNNNXXXX", X means random barcode, N means sequencing barcode')
    parser.add_argument('-i', '--input', required=True, help='input fastq file')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    return parser.parse_args()

def strip_barcode(pattern, seq):
    bc_len = len(pattern)
    strip_seq = seq[bc_len:]
    barcode = seq[:bc_len]
    return barcode, strip_seq


if __name__ == "__main__":   
    args = parse_args()
    bc_pattern = args.pattern
    prefix = args.input.split('fastq')[0]
    
    with open(args.input) as input_file, open(args.output,'w') as output_file, open(prefix+'bc', 'w') as bc_file:
        
        for tag, seq, strand, quality, in itertools.zip_longest(*[input_file]*4):
            barcode, seq = strip_barcode(args.pattern, seq.rstrip())
            tag = tag.rstrip() + '_' + barcode
            strand = strand.rstrip()
            quality = quality.rstrip()[len(barcode):]
            output_file.write(tag+'\n'+seq+'\n'+strand+'\n'+quality+'\n')
            bc_file.write(tag+'\t'+barcode+'\n')
    
        
            

