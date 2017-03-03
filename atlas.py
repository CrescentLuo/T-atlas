#!/usr/bin/env python

import time, os, glob, re, math, sys
import argparse
#APP_DIR=os.path.split(os.path.abspath(__file__))[0]
from parse import getFiles
import pysam, math
sys.path.insert(0, '/picb/extprog/biopipeline/data/database/20110414/UCSC/hg19')
from refFlat import RefFlat
debug = 0
# /picb/rnomics1/sszhu/singleWork/dataFromLi/wig_int_wald_rank_12_final.pl
#  $se_hat = sqrt( $p1_hat*(1-$p1_hat)/$n1 + $p2_hat*(1-$p2_hat)/$n2 );
#  $wald_stat = ($p1_hat - $p2_hat) / $se_hat;
# /picb/rnomics1/xiaoou/program/bamTobpkm.py
def pvalue(x):
    '''Cumulative distribution function for the standard normal distribution
    python
    import math
    x=2.33
    x=1.96
    1-(1.0 + math.erf(x / math.sqrt(2.0))) / 2.0
    x=-2.33
    (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0
    '''
    if x < 0:
        return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0
    else:
        return 1-(1.0 + math.erf(x / math.sqrt(2.0))) / 2.0
    # http://www.cyclismo.org/tutorial/R/pValues.html
    #R: 2*(1-pnorm(xbar,mean=a,sd=s/sqrt(20)))
    #2*pt(-abs(t),df=n-1)

def parallel(d, name):
    time.sleep(1)
    lock.acquire()
    print 'hello--'+str(name)
    print os.getppid(), '-----------', os.getpid()
    d[name] = name * name
    print d
    lock.release()

def getWaldScore(n1, N1, n2, N2):
    r1 = n1 / N1
    r2 = n2 / N2
    var = r1 * (1 - r1) / N1 + r2 * (1 - r2) / N2
    if var <= 0:
        return 0
    else:
        return (r2-r1)/math.sqrt(var)
def get_tissue_specific_score(vals):
### from paper 2014_Nature_Necsulea_The evolution of lncRNA repertoires and expression patterns in tetrapods;
    ma = max(vals)
    if ma <= 0: return "NA"
    nn = len(vals)
    if nn < 2: return "NA"
    score = 0
    for i in vals:
        score += (1 - (i / ma)) / (nn - 1)  ### Score(X)=sum((1-(x/max(X)))/n-1)
    return score

class myopen:
    def __init__(self, file):
        self.ff = open(file, 'w')
    def write(self, val):
        self.ff.write(val)
        if debug: print val,

def getCmps(cmps,flag2id):
    cmps_ids = []
    for cmp in re.split('\s*,\s*', cmps):
        a, b = cmp.split('/')
        assert a in flag2id and b in flag2id, "Samfiles in %s is not exists, Please check again."%cmp
        cmps_ids.append((flag2id[b], flag2id[a]))
    return cmps_ids

def summaryFiles(options):
    ### get refseqs
    refFlat=RefFlat(options.refFlat,symbols=options.genes)
    refFlat.getSym2ref()
    #return
    options.flags = []
    flag2id = {}
    ### get sam and size
    sam2size = {}
    sam2rlength = {}
    sams = []
    i = 0
    for samfile in options.infile:
        sam = pysam.AlignmentFile(samfile, "rb")
        tag = re.sub('\.[^.]*$', '', os.path.basename(samfile))
        #tag=re.sub('_unique_accepted_hits_sorted','',tag)
        options.flags.append(tag)
        flag2id[tag] = i
        i += 1
        ### get sam2size
        size = 0
        for line in os.popen("samtools idxstats %s"%samfile).readlines():
            size += int(line.split()[2])
        sam2size[sam] = size
        sams.append(sam)
        ### get sam2rlength
        for read in pysam.Samfile("%s"%samfile, "rb" ).fetch():
            sam2rlength[sam] = read.query_length
            break
        print "%s\tLength:%s"%(sam.filename, sam2rlength[sam])
    if not options.cmps:
        if len(sams) > 1:
            options.cmps = [(0, 1)]
    else:
        options.cmps = getCmps(options.cmps, flag2id)

    ff = myopen(options.outfile)
    ff.write("Symbol\tLoc\tStrand\tLength\tNumberOfExons\tTheNthExons\t"+"\t".join(options.flags)+"\tSpecificity_Score\t")
    if options.cmps:
        ff.write("\t".join(["WaldScore(%s/%s)\tP_value(%s/%s)\tFold((Y+1)/(X+1))"%(options.flags[b],options.flags[a], options.flags[b],options.flags[a]) for a, b in options.cmps]))
    ff.write("\n")
    N = len(refFlat.sym2ref.values())
    ii = 0
    for ref in refFlat.sym2ref.values():
        ii += 1
        #if ii<20890: continue
        if not ii % 1000: print "%s/%s processed"%(ii, N)
        ref.getRefSeq2eles()
        #if len(ref.eles.exons)<5: continue
        # print ref.sym,ref.data,ref.eles.exons
        sam2values = {}
        sam2lens = {}
        sam2HPBs = {}
        starts = []
        ends = []
        for aa, bb in ref.eles.exons:
            aa, bb = int(aa), int(bb)
            starts.append(aa)
            ends.append(bb)
            values=[]
            for sam in sams:
                n=0
                for aread in sam.fetch(ref.chr,aa,bb):
                    n+=aread.get_overlap(aa,bb)
                values.append(n)
                sam2values.setdefault(sam,[]).append(n)
                sam2lens.setdefault(sam,[]).append(bb-aa)
                sam2HPBs.setdefault(sam,[]).append(n*1.0*10**9/(sam2rlength[sam]*sam2size[sam]*(bb-aa)))
        ### get sam2sum and sam2len
        sam2sum={}
        sam2len={}
        for sam in sams:
            sam2sum[sam]=sum(sam2values[sam])
            sam2len[sam]=sum(sam2lens[sam])
        ### get sam2HPB and HPBs
        sam2HPB={}
        HPBs=[]
        for sam in sams:
            HPB=sam2sum[sam]*1.0*10**9/(sam2size[sam]*sam2rlength[sam]*sam2len[sam])
            sam2HPB[sam]=HPB
            HPBs.append(HPB)
        if sum(HPBs)<=1: continue
        ### print   the total exon marked as 0 exon
        ff.write("%s\t%s\t%s\t%s\t%s\t0\t%s\t%s\t"%(ref.sym,"%s:%s-%s"%(ref.chr,ref.la,ref.lb),ref.sign,ref.lb-ref.la,len(ref.eles.exons),"\t".join(map(str,HPBs)),get_tissue_specific_score(HPBs)))
        valds=[]
        if options.cmps:
            for ia,ib in options.cmps:
                sa=sams[ia]
                sb=sams[ib]
                vald_t=getWaldScore(sam2sum[sa]*1.0/sam2rlength[sa],sam2size[sa],sam2sum[sb]*1.0/sam2rlength[sb],sam2size[sb])
                valds.extend([vald_t,pvalue(abs(vald_t)),"%s"%((sam2HPB[sb]+1)*1.0/(sam2HPB[sa]+1))])
        ff.write("\t".join(map(str,valds))+"\n")
        ### print each exons
        ### p_values = scipy.stats.norm.sf(z_scores) #one-sided
        ### p_values = scipy.stats.norm.sf(z_scores)*2 #twosided
        ### >>>normpdf(7,5,5)  
        ### 0.073654028688865794
        ### >>> norm(5,5).pdf(7)
        ### >>> 0.073654028060664664
        nn=len(sam2values[sams[0]])
        if nn<=1: continue
        for i in range(nn):
            vals=[sam2HPBs[sam][i] for sam in sams]
            ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(ref.sym,"%s:%s-%s"%(ref.chr,starts[i],ends[i]),ref.sign,ends[i]-starts[i],nn,i+1,"\t".join(map(str,vals)),get_tissue_specific_score(vals)))
            valds=[]
            if options.cmps:
                for ia,ib in options.cmps:
                    sa=sams[ia]
                    sb=sams[ib]
                    vald=getWaldScore(sam2values[sa][i]*1.0/sam2rlength[sa],sam2size[sa],sam2values[sb][i]*1.0/sam2rlength[sb],sam2size[sb])
                    valds.extend([vald,pvalue(abs(vald)),"%s"%((sam2HPBs[sb][i]+1)*1.0/(sam2HPBs[sa][i]+1))])
            ff.write("\t".join(map(str,valds))+"\n")
        #break

def rprint(f,val):
    f.write(val)
    print val,


if __name__ == '__main__':
    '''
    parseBam.py -i '*.bam' -o tt.txt
    parseBam.py -i "*.bam" -c "HCT116_CCAT1_kd/HCT116_scramble" -o tt.txt
    parseBam.py -h
    '''
    #python
    #import argparse
    des='''
This program is used to calculate BPKM values for each merged gene in refFlat file and different comparisons.
Below is the perl script to transform knowGene.txt to refFlat format:
perl -F\'\\t\' -ane \'if($#F<=9){$F[4]=~s/ /_/g;$a{$F[0]}=$F[4]}else{$"="\\t";print "$a{$F[0]}\\t@F[0..9]\\n"}\' kgXref.txt knownGene.txt > knownRef.txt
\nWald test score: r1=n1/N1; r2=n2/N2; var=r1*(1-r1)/N1+r2*(1-r2)/N2; score=(r2-r1)/math.sqrt(var)
Specificity score: sum((1-(x/max(X)))/n-1)
    '''
    parser = argparse.ArgumentParser(description=des,
            epilog='eg: %(prog)s -i "*.bam" -c "HCT116_CCAT1_kd/HCT116_scramble,PA1_SL5AC_kd/PA1_scramble_1" -o tt.txt',
            prefix_chars="-+", fromfile_prefix_chars="@",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', dest='infile',  help='input file name. "," seprated bam files, support "*" and "|", such as -i "*.bam|tt/*.bam,mm/*.bam"', required =True)
    #parser.add_argument('-f', '--fasta', dest='fasta',  help='genome fasta file.',default='/picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/seqs/hg19.fa')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name, default="output.txt"', default='output.txt')
    refs=['"-r /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/refFlat.txt"',
        "Other useful annotations:",
        "-r /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/knownRef.txt",
        "-r /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/20131128/refFlat.txt",
        "-r /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/20131128/knownRef.txt"]
    parser.add_argument('-r', '--refFlat', dest='refFlat', help='gene annotation file from UCSC, default='+"\n".join(refs), default='/picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/refFlat.txt')
    parser.add_argument('-g', '--genes', dest='genes',  help='Only look at the selected genes, "," seprated comparison, default=None')
    parser.add_argument('-c', '--cmps', dest='cmps',  help='Compare between two samples, "," seprated comparison, such as -c "HCT116_CCAT1_kd/HCT116_scramble,PA1_SL5AC_kd/PA1_scramble_1", default=file[1]/file[0]')
    options= parser.parse_args()
    print '###Parameters:'
    for key,val in options.__dict__.items():
        print '%s:\t%s'%(key,val)
    print '###Parameters:\n',
    options.infile=getFiles(options.infile)
    if options.genes:
        options.genes=re.split("\s*, \s*", options.genes)
    summaryFiles(options)

   