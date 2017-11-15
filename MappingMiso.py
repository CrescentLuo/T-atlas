import pandas as pd
import argparse
from BCBio import GFF

parser = argparse.ArgumentParser(description='This script is written for mapping MISO events id to gene annotation')
parser.add_argument('-i', '--input', help="Input from miso filter output")
parser.add_argument('-g', '--gff', help='gene annotaion gff file')
args = parser.parse_args()

limits = dict(gff_type = ['gene','mRNA'])
gff_handle = open(args.gff)
gff_dict = dict()
for rec in GFF.parse(gff_handle, limit_info=limits):
    for gene_feature in rec.features:
        print gene_feature
        print gene_feature.id, gene_feature.qualifiers['gsymbol'][0]
        gff_dict[gene_feature.id] = [gene_feature.qualifiers['ensg_id'][0],gene_feature.qualifiers['gsymbol'][0]]

miso_bf = pd.read_table(args.input)
for index,event in miso_bf.iterrows():
    print event['event_name'], gff_dict[event['event_name']][1], event['diff']

