import sys
import argparse
import pysam
import pybedtools
from pybedtools import BedTool
import re
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract exons from bed file")
    parser.add_argument('--anno', required=True,
                        help='input annotation in gtf format')
    parser.add_argument('--gene_type', default=None, help='specific genetype to exact')
    parser.add_argument('--mid_exon_only', action='store_true', help="only extract middle exons" )
    parser.add_argument('-o', '--output', default='./exons.bed', help="output prefix for files")
    args = parser.parse_args()
    return args


def get_transcript_id(transcript):
  ts_id = transcript.attrs['transcript_id']
  return ts_id

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


def get_exon(anno_file, gene_type, mid=False):
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
                 if get_transcript_id(e) == ts_id]
        exons.sort()

        if not re.match("^chr((\d+)|X|Y)$",ts.chrom):
            continue

        name_format_str = "{gene_id}_{gene_name}"

        for i, exon in enumerate(exons):
          if mid and (i==0 or i==len(exons)-1):
            continue
          else:
            yield pybedtools.Interval(
                ts.chrom, exon.start, exon.end,
                name_format_str.format(**get_gene_id_and_name(ts)),
                "exon", ts.strand
            )


if __name__ == "__main__":
    args = parse_args()
    with open(BedTool._tmp(), "w") as fh:
      for i, exon in enumerate(get_exon(args.anno, args.gene_type, args.mid_exon_only)):
        fh.write(str(exon))
      exons = BedTool(fh.name)
      exons.sort().saveas(args.output)
      srt_output_fn = ".".join(args.output.split('.')[:-1]) + ".srt.bed" 
      with open(srt_output_fn, "w") as srt_output:
        subprocess.call(["sort", "-k","1,1","-k", "2,2n", "-k", "3,3n", "-k", "6,6","-u", args.output], stdout=srt_output_fn)
