import sys
import argparse
import pysam
import pybedtools
from pybedtools import BedTool
import tqdm
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract introns from bed file")
    parser.add_argument('--anno', required=True,
                        help='input annotation in gtf format')
    parser.add_argument('--gene_type', default=None,
                        help='specific genetype to exact')
    parser.add_argument(
        '-o', '--output', default='./introns.bed', help="output prefix for files")
    args = parser.parse_args()
    return args


def get_transcript_id(transcript):
  ts_id = transcript.attrs['transcript_id']
  return ts_id


def get_interval_iter(anno, feature, gene_type=None):
  print("searching for {}".format(feature), file=sys.stderr)
  if gene_type == None:
    return BedTool((interval for interval in anno if interval[2] == feature)).saveas()
  else:
    return BedTool((interval for interval in anno if interval[2] == feature and interval.attrs['gene_type'] == gene_type)).saveas()


def get_gene_id_and_name(transcript):
  gene_id = transcript.attrs['gene_id']
  gene_name = transcript.attrs['gene_name']
  return dict(gene_id=gene_id, gene_name=gene_name)


def get_intron(anno_file, gene_type):
    anno = pybedtools.BedTool(anno_file)
    print("Read input annotaion file", file=sys.stderr)
    anno_srt = anno.sort()
    print("Sort input annotation file", file=sys.stderr)
    exon_iter = get_interval_iter(anno_srt, "exon", gene_type)
    trnascript_iter = get_interval_iter(anno_srt, "transcript", gene_type)
    print("Fetch iters of exons and transcripts", file=sys.stderr)
    exon_intervals = pybedtools.IntervalFile(exon_iter.fn)
    print("Generate intervals", file=sys.stderr)
    for ts in trnascript_iter:
        if not re.match("^chr((\d+)|X|Y)$",ts.chrom):
          continue
        ts_id = get_transcript_id(ts)
        exons = [e for e in exon_intervals.search(
            ts, same_strand=True) if get_transcript_id(e) == ts_id]
        exons.sort()
        

        name_format_str = "{gene_id}_{gene_name}"

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
