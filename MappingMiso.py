import argparse
from BCBio import GFF
parser = argparse.ArgumentParser(description='This script is written for mapping MISO events id to gene annotation')
parser.add_argument('-i', '--input', help="Input from miso filter output")
parser.add_argument('-g', '--gff', help='gene annotaion gff file')
args = parser.parse_args()

limits = dict(gff_type = ['gene','mRNA'])
gff_handle = open(args.gff)
for rec in GFF.parse(gff_handle, limit_info=limits):
    for gene_feature in rec.features:
        print gene_feature

