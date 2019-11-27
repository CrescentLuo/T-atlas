import sys
import argparse
import pysam
import pybedtools
from pybedtools import BedTool
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract exons and introns from bed file")
    parser.add_argument('--anno', required=True,
                        help='input annotation in gtf format')
    parser.add_argument('--gene-type', default=None, help='specific genetype to exact')
    parser.add_argument('-o', '--output', default='./introns.bed', help="output prefix for files")
    args = parser.parse_args()
    return args


def get_transcript_id(transcript):
  if 'transcript_id' in transcript.attrs:
    ts_id = transcript.attrs['transcript_id']
  else:
    ts_id = re.match(r'transcript:(.*)', transcript.attrs['ID']).group(1)
  return ts_id


def get_parent_transcript_id(exon):
  if 'transcript_id' in exon.attrs:
    return exon['transcript_id']
  else:
    return re.match(r'transcript:(.*)', exon.attrs['Parent']).group(1)


def get_interval_iter(anno, feature, gene_type=None):
  if gene_type == None:
    return BedTool((interval for interval in anno if interval[2] == feature)).saveas()
  else:
    return BedTool((interval for interval in anno if interval[2] == feature and interval.attrs['gene_type'] == gene_type)).saveas()


def get_gene_id_and_name(transcript):
  if 'gene_id' in transcript.attrs:
    gene_id = transcript.attrs['gene_id']
  else:
    gene_id = re.match(r'gene:(.*)', transcript.attrs['Parent']).group(1)

  if 'gene_name' in transcript.attrs:
    gene_name = transcript.attrs['gene_name']
  elif 'Name' in transcript.attrs:
    gene_name = re.match(r'(.*)-(\d+)', transcript.attrs['Name']).group(1)
  else:
    gene_name = ""
  return dict(gene_id=gene_id, gene_name=gene_name)


def get_intron(anno_file, gene_type):
    anno = pybedtools.BedTool(anno_file)
    anno_srt = anno.sort()

    ts_dict = dict()
    ts_exon = dict()
    exon_iter = get_interval_iter(anno_srt, "exon", gene_type)
    trnascript_iter = get_interval_iter(anno_srt, "transcript", gene_type)

    exon_intervals = pybedtools.IntervalFile(exon_iter.fn)

    for ts in trnascript_iter:
        ts_id = get_transcript_id(ts)

        exons = [e for e in exon_intervals.search(ts, same_strand=True)
                 if get_parent_transcript_id(e) == ts_id]
        exons.sort()

        if not ts.chrom.startswith('chr'):
            ts.chrom = "chr{}".format(ts.chrom)

        name_format_str = "{gene_id} ({gene_name})"

        for i, exon in enumerate(exons):
          if i < len(exons) - 1:
            istart = exon.end
            iend = exons[i + 1].start
            assert istart < iend
            yield pybedtools.Interval(
                ts.chrom, istart, iend,
                name_format_str.format(**get_gene_id_and_name(ts)),
                "intron", ts.strand
            )


if __name__ == "__main__":
    args = parse_args()
    with open(BedTool._tmp(), "w") as fh:
      for i, intron in enumerate(get_intron(args.anno, args.gene_type)):
        fh.write(str(intron))
      introns = BedTool(fh.name)
      introns.sort().merge(s=True, c=(4, 5, 6), o='distinct').saveas(args.output)
