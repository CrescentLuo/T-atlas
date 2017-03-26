import re, sys, os
import warnings
import pysam
from genomeInfor import genome2seq
from Bio.Seq import Seq
# fasta=genome2seq().get_genome_seq("hg19")
def custom_formatwarning(msg, *a):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning
#warnings.warn("achtung")

dataDir = os.path.split(__file__)[0]

class ref2infor(dict):
    def __init__(self, file=os.path.join(dataDir, 'refFlat.txt')):
        ref2infor = {}
        for line in open(file):
            xl = line.rstrip().split('\t')
            if len(xl) < 3: continue
            if not re.match('chr‘\w’{,2}$',xl[2]): continue
            #### use loc as unique id
            ref=Ref2eles(xl)
            ref2infor[xl[1]]=ref
        dict.__init__(self,ref2infor)

class Refs(list):
    def __init__(self,file=os.path.join(dataDir,'refFlat.txt')):
        refs=[]
        for line in open(file):
            xl=line.rstrip().split('\t')
            if len(xl)<3: continue
            if not re.match('chr\w{,2}$',xl[2]): continue
            #### use loc as unique id
            ref=Ref2eles(xl)
            refs.append(ref)
        list.__init__(self,refs)
    #### define iterator
    #def __iter__(self):
    #    self.n=-1
    #    return self
    #def next(self):
    #    self.n+=1
    #    if self.n>=len(self.refs): raise StopIteration
    #    return self.refs[self.n]

    #def __repr__(self,refs=['NM_175061']):
    #    frefs=[ref for ref in refs if ref in self]
    #    if len(frefs)<1: return ''
    #    for ref in frefs:
    #        return str(ref)
class Ref2eles:   #### parse elements of each refseq
    def __init__(self,xl):
        self.data = xl
        self.ref = xl[1]
        self.ele2locs = {}
        self.exonsAdn3k = ''
        self.exons = ''
        self.introns = ''    ### updated by sszhu1007@gmail.com 11.11.17 12:11:20 ###
        self.dn3k = ''
        self.promoter = ''
        self.chr = xl[2]
        self.sym = xl[0]
        self.sign = xl[3]
        self.la = xl[4]
        self.lb = xl[5]
        self.TSS, self.TTS = xl[4], xl[5]
        self.ele2locs['L_UTR'] = [(xl[4], xl[6])]
        self.ele2locs['R_UTR'] = [(xl[7], xl[5])]
        self.ele2locs['5_UTR'], self.ele2locs['3_UTR'] = \
            self.ele2locs['L_UTR'], self.ele2locs['R_UTR']
        self.ele2locs['promoter'] = [(str(int(xl[4])-2000), xl[4])]
        self.promoter = [(str(int(xl[4]) - 2000), str(int(xl[4]) + 2000))]
        self.ele2locs['promoter_body'] = [(str(int(xl[4]) - 2000), xl[5])]
        self.dn3k = (xl[5], str(int(xl[5]) + 3000))
        if xl[3] == '-':
            self.TSS,self.TTS = xl[5],xl[4]
            self.ele2locs['5_UTR'], self.ele2locs['3_UTR'] = \
                self.ele2locs['R_UTR'],self.ele2locs['L_UTR']
            self.ele2locs['promoter']=[(xl[5],str(int(xl[5])+2000))]
            self.promoter=[(str(int(xl[5])-2000),str(int(xl[5])+2000))]
            self.ele2locs['promoter_body']=[(xl[4],str(int(xl[5])+2000))]
            st=int(xl[4])-3000
            if st<0: st=0
            self.dn3k=(str(st),xl[4])
        las=xl[9].split(',')[:-1]
        lbs=xl[10].split(',')[:-1]
        self.exons=zip(las,lbs)
        self.ele2locs['exons']=zip(las,lbs)
        self.introns=zip(lbs,las[1:])    ### updated by sszhu1007@gmail.com 11.11.17 12:18:22 ###
        self.ele2locs['introns']=zip(lbs,las[1:])
        if xl[3]=='-':
            self.exons.reverse()
            self.ele2locs['exons'].reverse()
            self.introns.reverse()    ### updated by sszhu1007@gmail.com 11.11.17 12:19:20 ###
            self.ele2locs['introns'].reverse()
            self.exonsAdn3k=[self.dn3k]+self.exons
        else:
            self.exonsAdn3k=self.exons+[self.dn3k]
    def __eq__(self):
        return self.ref
    def getEles(self,tag,locs):
        aa=int(self.__dict__[tag])
        olocs=[]
        for a,b in locs:
            if self.sign=='+':
                olocs.append((aa+int(a),aa+int(b)))
            else:
                olocs.append((aa-int(b),aa-int(a)))
        return olocs
    def getRegionEles(self,tag1,tag2,locs):
        aa=int(self.__dict__[tag1])
        bb=int(self.__dict__[tag2])
        olocs=[]
        for a,b in locs:
            if self.sign=='+':
                olocs.append((aa+int(a),bb+int(b)))
            else:
                olocs.append((bb-int(b),aa-int(a)))
        return olocs
    def __repr__(self):
        ele2locs=self.ele2locs
        chr=self.chr
        sign=self.sign
        ref=self.ref
        out=[]
        for ele,locs in ele2locs.items():
            n=len(locs)
            if n<2:
                out.append('\t'.join([chr]+list(locs[0]))+'\t%s_%s_%s\t%s'%(self.sym,ref,ele,sign))
            else:
                for id in range(n):
                    out.append('\t'.join([chr]+list(locs[id]))+'\t%s_%s_%s_%s\t%s'%(self.sym,ref,id+1,ele,sign))
        return '\n'.join(out)

class GeneInfor:
    def __init__(self,file='',taxId=9606):
        if not (os.path.isfile(file) and os.path.getsize(file) > 100): file="/picb/extprog/biopipeline/data/database/GO20130411/gene_info_%s"%taxId
        if not (os.path.isfile(file) and os.path.getsize(file) > 100): file="/picb/extprog/biopipeline/data/database/20110414/NCBI_GO/gene_info_%s"%taxId
        if not (os.path.isfile(file) and os.path.getsize(file) > 100): 
            print "No gene info file: %s"%file
            return
        self.refFile=file
        self.taxId=taxId
        self.sym2Synonyms={}
        self.sym2Description={}
        self.geneid2sym={}
        self.geneid2Synonyms={}
        self.geneid2Description={}
        for line in open(file):
            '''
            #Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status Other_designations Modification_date (tab is used as a separator, pound sign - start of a comment)
            9606	1	A1BG	-	A1B|ABG|DKFZp686F0970|GAB|HYST2477	HGNC:5|MIM:138670|Ensembl:ENSG00000121410|HPRD:00726	19	19q13.4	alpha-1-B glycoprotein	protein-coding	A1BG	alpha-1-B glycoprotein	O	alpha-1B-glycoprotein	20110403
            9606	2	A2M	FWP007	CPAMD5|DKFZp779B086|S863-7	HGNC:7|MIM:103950|Ensembl:ENSG00000175899|HPRD:00072	12	12p13.31	alpha-2-macroglobulin	protein-coding	A2M	alpha-2-macroglobulin	O	C3 and PZP-like alpha-2-macroglobulin domain-containing protein 5|alpha-2-M	20110403
            9606	3	A2MP1	-	A2MP	HGNC:8	12	12p13.3-p12.3	alpha-2-macroglobulin pseudogene 1	pseudo	A2MP1	alpha-2-macroglobulin pseudogene 1	O	-	20101101
            '''
            xl=line.rstrip().split('\t')
            if len(xl)<3: continue
            self.sym2Synonyms[xl[2]]=xl[4]
            self.sym2Description[xl[2]]=xl[8]
            self.geneid2sym[xl[1]]=xl[2]
            self.geneid2Synonyms[xl[1]]=xl[4]
            self.geneid2Description[xl[1]]=xl[8]
    def getGeneInfor(self,tags,gene):
        out=[]
        for tag in tags:
            out.append(self.__dict__.get(tag,"").get(gene,""))
        return out

class RefFlat:
    def __init__(self,file=os.path.join(dataDir,'refFlat.txt'),taxId=0,symbols=''):
        #self.locs=set()
        self.chrs=set()
        self.chr2locs={}
        self.refSeqs=set()
        self.ref2info={}    ### updated by sszhu1007@gmail.com 11.10.11 08:31:21 ###
        self.chr2bound=0
        self.refFile=file
        self.sym2refs={}
        self.ref2chrs={}
        for line in open(file):
            '''
            SNORA1	NR_003026	chr11	-	93465169	93465299	93465299	93465299	1	93465169,	93465299,
            SNORA59A	NR_003025	chr1	+	12567299	12567451	12567451	12567451	1	12567299,	12567451,
            DDX3Y	NM_004660	chrY	+	15016698	15032390	15016847	15030034	17	15016698,15019447,15021270,15023750,15024638,15024874,15025629,15026475,15026795,15026978,15027541,15027794,15028172,15028428,15028818,15029314,15029954,	15016892,15019505,15021318,15023880,15024794,15024974,15025765,15026561,15026894,15027139,15027686,15027939,15028354,15028546,15028972,15029454,15032390,
            '''
            xl=self.getFormatxl(line,file)
            #xl=line.rstrip("\n").split('\t')
            if symbols and (xl[0] not in symbols): continue
            if len(xl)<3: continue
            if not re.match('chr\w{,2}$', xl[2]): continue
            #### use loc as unique id
            #loc='%s:%s-%s:%s'%(xl[2],xl[4],xl[5],xl[3])
            #if loc in self.locs: continue
            #loc='%s:%s:%s-%s:%s'%(xl[1],xl[2],xl[4],xl[5],xl[3])    ### updated by sszhu1007@gmail.com 12.09.25 12:31:01 ###
            self.chrs.add(xl[2])
            self.ref2chrs.setdefault(xl[1],set()).add(xl[2])         ### remove isoforms in different chromosomes
            if len(self.ref2chrs[xl[1]])>1: continue
            #self.locs.add(loc)
            obj=RefSeq(xl)
            self.chr2locs.setdefault(xl[2],set()).add(obj)
            self.sym2refs.setdefault(obj.sym,[]).append(obj)
            #self.chr2locs.setdefault(xl[2],[]).append(obj)
            #self.allLocs.add(obj)    ### updated by sszhu1007@gmail.com 11.11.17 11:23:20 ###
            self.refSeqs.add(obj)
            self.ref2info[xl[1]]=obj
        ### delete without geneid in allLocs and get geneid2loc
        if taxId: self.updateGeneInfor(taxId)

    def getFormatxl(self,line,file):
        xl=line.rstrip("\n").split('\t')
        n=len(xl)
        if re.search('.txt$',file):     #### for refFlat.txt
            return xl
            #SNORA1	NR_003026	chr11	-	93465169	93465299	93465299	93465299	1	93465169,	93465299,
        elif re.search('.bed$',file):   #### for bed files
            tag='%s:%s-%s' % (xl[0], xl[1], xl[2])
            if n==3 and re.search("^chr",xl[0]):   #### 3 column bed
                ### name name chr sign la lb lb lb e_number e_sts e_ends
                return [tag,tag,xl[0],"+",xl[1],xl[2],xl[2],xl[2],1,xl[1]+",",xl[2]+","]
            elif n==6 and re.search("^chr",xl[0]):   #### 6 column bed
                # chr9    136325086       136335909       CACFD1  NM_001242370    +
                # chrX    51804922        51812368        MAGED4B NM_001242362    -
                return [xl[3]+"_"+xl[4],tag,xl[0],xl[5],xl[1],xl[2],xl[2],xl[2],1,xl[1]+",",xl[2]+","]
            elif n==12 and re.search("^chr",xl[0]):   #### 12 column bed
                size = xl[10].split(',')
                offset = xl[11].split(',')
                starts = [int(xl[1]) + int(i) for i in offset]
                start =",".join(map(str,starts))+","
                ends =[starts[i] + size[i] for i in xrange(len(offset))]
                end =",".join(map(str,ends))+","
                return [xl[3],tag,xl[0],xl[5],xl[1],xl[2],xl[6],xl[7],xl[9],start,end]
            else:
                sys.exit('Error: not suported bed format')
                # /picb/extprog/biopipeline/bin/splicing.pl  ### for junction reads
                #my $chr=$F[0];
                #my @sizes=split /,/,$F[10];
                #my @starts=split /,/,$F[11];
                #my ($la,$lb)=($F[1]+$starts[0]+$sizes[0],$F[1]+$starts[1]);
                ##chr1	12214	12681	JUNC00000001	1	+	12214	12681	255,0,0	2	13,87	0,380
        else:
            sys.exit('Error: not suported format!!')

    def getSym2ref(self):
        self.sym2ref={}
        warns=[]
        for sym in self.sym2refs:
            nn=len(self.sym2refs[sym])
            ##### no merge for singleton ref
            if nn==1:
                self.sym2ref[sym]=self.sym2refs[sym][0]
                continue
            ##### merge multi refs for each gene
            ref=self.sym2refs[sym][0]
            la,lb,lc,ld=ref.la,ref.lb,ref.lc,ref.ld
            starts=[]
            ends=[]
            starts.extend(ref.starts)
            ends.extend(ref.ends)
            chrs=set()
            chrs.add(ref.chr)
            ## merge la,lb,lc,ld
            for ref in self.sym2refs[sym][1:]:
                if ref.la< la: la=ref.la
                if ref.lc< lc: lc=ref.lc
                if ref.lb> lb: lb=ref.lb
                if ref.ld> ld: ld=ref.ld
                starts.extend(ref.starts)
                ends.extend(ref.ends)
                chrs.add(ref.chr)
            if len(chrs)>1:
                for ref in self.sym2refs[sym]:
                    warns.append("\t".join(map(str,ref.data)))
                continue
            ## merge starts and ends for all exons
            starts,ends=self.mergeExons(starts,ends)
            xl=map(str,[ref.sym,"merge_%s"%nn,ref.chr,ref.sign,la,lb,lc,ld,len(starts),",".join(map(str,starts)),",".join(map(str,ends))])  #### add "," in mergeExons
            ## construct new RefSeq object
            self.sym2ref[sym]=RefSeq(xl)
        n=len(warns)
        if n>0:
            genes=set()
            for data in warns:
                xl=data.rstrip().split('\t')
                genes.add(xl[0])
            nn=len(genes)
            warnings.warn("We skiped %s genes, which located on multiple chromosomes!!!!\n"%nn+"\n".join(warns))

    def printSym2merge(self,out="sym2merge.txt"):
        self.getSym2ref()
        fout=open(out,"w")
        for sym in self.sym2ref:
            fout.write("\t".join(self.sym2ref[sym].data)+"\n")

    def mergeExons(self,starts,ends):
        #### get exons with type
        exons=[]     # position of [ or ] boundary  #  -1 => [   +1 => ]   
        exons.extend(zip(starts,[-1]*len(starts)))  #### [
        exons.extend(zip(ends,[1]*len(ends)))       #### ]
        exons.sort()
        #### get interval
        pre_bracket=0        #bracket position
        pre_type=0           #[ => -1    ] => +1
        count=0
        starts=[]
        ends=[]
        for bracket,type in exons:
            if bracket != pre_bracket and count !=0:
                #print "%s\t%s\t%s"%(pre_bracket,bracket,count)
                starts.append(pre_bracket)
                ends.append(bracket)
            count=count-type
            pre_bracket=bracket
            pre_type=type
        #if len(starts)<=1:
        starts.append('')                        #### to simulate refFlat.txt
        ends.append('')
        return starts,ends

#    def updateGeneInfor(self,taxId):
#        ### delete without geneid
#        dels=set()
#        import sys
#        sys.path.extend(['/home/sszhu/sszhu10/projects/gene_data_processing','/home/sszhu/sszhu10/projects/gene_data_processing/hg18_20100305'])
#        import geneMapping
#        #taxId=9606
#        #taxId=10090
#        geneid2name=geneMapping.geneid2name(tax=taxId)
#        name2geneid=geneid2name.name2geneid
#        geneid2des=geneid2name.geneid2des
#        geneid2name=geneid2name.geneid2name
#        self.geneid2loc={}
#        for loc in self.allLocs:
#            if loc.sym in name2geneid:
#                #### get geneid from symbol
#                id=name2geneid[loc.sym]
#                loc.geneid=str(id)
#                loc.des=geneid2des[id]
#                if not id in self.geneid2loc or self.geneid2loc[id].geneLength<loc.geneLength:
#                    self.geneid2loc[id]=loc
#            else:
#                dels.add(loc)
#                #print loc
#        self.allLocs=self.allLocs - dels
    def getExonBoundary(self):
        self.chr2bound=1
        self.chr2boundary={}
        for line in open(self.refFile):    ##### use all refSeqs
            '''
            SNORA1	NR_003026	chr11	-	93465169	93465299	93465299	93465299	1	93465169,	93465299,
            SNORA59A	NR_003025	chr1	+	12567299	12567451	12567451	12567451	1	12567299,	12567451,
            DDX3Y	NM_004660	chrY	+	15016698	15032390	15016847	15030034	17	15016698,15019447,15021270,15023750,15024638,15024874,15025629,15026475,15026795,15026978,15027541,15027794,15028172,15028428,15028818,15029314,15029954,	15016892,15019505,15021318,15023880,15024794,15024974,15025765,15026561,15026894,15027139,15027686,15027939,15028354,15028546,15028972,15029454,15032390,
            '''
            xl=line.rstrip().split('\t')
            chr=xl[2]
            self.chr2boundary.setdefault(chr,set())
            starts=map(int,xl[9].split(',')[:-1])
            ends=[x-1 for x in map(int,xl[10].split(',')[:-1])]
            self.chr2boundary[chr].update(starts)
            self.chr2boundary[chr].update(ends)

    def nearestDistFromExonBoundary(self,loc):
        #### calculate distance from exon, consider strand -:up, +:down;
        ### get self.chr2boundary
        if not self.chr2bound:
            self.getExonBoundary()
        tloc=Loc(loc)
        ##ldis=rdis=sys.maxint
        dis=100000000
        loct='NA'
        if not tloc.chr in self.chrs:
            return 'NA','NA'
        for loc in self.chr2boundary[tloc.chr]:
            if tloc.la<=loc<tloc.lb:
                return 0,loc
            if tloc.la>loc: ###loc locate tloc's left
                ddis=tloc.la - loc
            else:              ### loc locate tloc's right
                ddis=loc - tloc.lb+1
            if abs(ddis)<abs(dis):
                dis=ddis
                loct=loc
        return dis,loct

    def nearestDistFromTSS(self,loc):
        #### calculate distance from TSS, consider strand -:up, +:down;
        tloc=Loc(loc)
        ##ldis=rdis=sys.maxint
        dis=100000000
        loct='na'
        if not tloc.chr in self.chrs:
            return 'NA'
        for loc in self.chr2locs[tloc.chr]:
            if tloc.la<=loc.TSS<=tloc.lb:
                return 0
            if tloc.la>loc.TSS: ###loc locate tloc's left
                ddis=tloc.la - loc.TSS
                if abs(ddis)<abs(dis):
                    dis=ddis
                    loct=loc
                    if loc.sign=='-':
                        dis=0-ddis
            else:              ### loc locate tloc's right
                ddis=tloc.lb - loc.TSS    ### checked by sszhu1007@gmail.com 12.03.08 10:02:39 ###
                if abs(ddis)<abs(dis):
                    dis=ddis
                    loct=loc
                    if loc.sign=='-':
                        dis=0-ddis
        return dis
        return dis,loct,loct.loc
    def nearbyLocs(self,loc):
        tloc=Loc(loc)
        ##ldis=rdis=sys.maxint
        ldis=rdis=100000000
        lloc=rloc='NA'
        overLocs=set()
        if not tloc.chr in self.chrs:
            return ldis,lloc,overLocs,rdis,rloc
        for loc in self.chr2locs[tloc.chr]:
            if loc.la<=tloc.la<=loc.lb or loc.la<=tloc.lb<=loc.lb \
                or tloc.la<=loc.la<=tloc.lb:
                overLocs.add(loc)
                continue
            if tloc.la>loc.lb: ###loc locate tloc's left
                lldis=tloc.la - loc.lb
                if lldis<ldis:
                    ldis=lldis
                    lloc=loc
            else:              ### loc locate tloc's right
                rrdis=loc.la - tloc.lb
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if rrdis<rdis:
                    rdis=rrdis
                    rloc=loc
        return ldis,lloc,overLocs,rdis,rloc

    def siteNearbyGenes(self,site):
        chr,mid=site.split(":")
        mid=int(mid)
        ##ldis=rdis=sys.maxint
        ldis=rdis=100000000
        lloc=rloc='NA'
        overLocs=set()
        if not chr in self.chrs:
            return ldis,lloc,overLocs,rdis,rloc
        for loc in self.chr2locs[chr]:
            if loc.la<=mid<=loc.lb:
                overLocs.add(loc)
                continue
            if mid>loc.lb: ###loc locate tloc's left
                lldis=mid - loc.lb
                if lldis<ldis:
                    ldis=lldis
                    lloc=loc
            else:              ### loc locate tloc's right
                rrdis=loc.la - mid
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if rrdis<rdis:
                    rdis=rrdis
                    rloc=loc
        return ldis,lloc,overLocs,rdis,rloc
    def nearbyTSS(self,loc):
        tloc=Loc(loc)
        ##ldis=rdis=sys.maxint
        ldis=rdis=100000000
        lloc=rloc='NA'
        overLocs=set()
        if not tloc.chr in self.chrs:
            return ldis,lloc,overLocs,rdis,rloc
        for loc in self.chr2locs[tloc.chr]:
            if tloc.la<=loc.TSS<=tloc.lb:
                overLocs.add(loc)
                continue
            if tloc.la>loc.TSS: ###loc locate tloc's left
                lldis=tloc.la - loc.TSS
                if lldis<ldis:
                    ldis=lldis
                    lloc=loc
            else:              ### loc locate tloc's right
                rrdis=loc.TSS - tloc.lb
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if rrdis<rdis:
                    rdis=rrdis
                    rloc=loc
        if len(overLocs)<1: overLocs=''
        return ldis,lloc,overLocs,rdis,rloc
    def nearTSS(self,loc,up=1000,dn=1000):
        tloc=Loc(loc)
        overLocs=set()
        if not tloc.chr in self.chrs:
            return overLocs
        for loc in self.chr2locs[tloc.chr]:
            if tloc.la<=loc.TSS<=tloc.lb:
                overLocs.add(loc)
                continue
            if tloc.la>loc.TSS: ###loc locate tloc's left
                lldis=tloc.la - loc.TSS
                if loc.sign=='+':
                    ldis=dn
                else:
                    ldis=up
                if lldis<=ldis:
                    overLocs.add(loc)
            else:              ### loc locate tloc's right
                rrdis=loc.TSS - tloc.lb
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if loc.sign=='+':
                    rdis=up
                else:
                    rdis=dn
                if rrdis<=rdis:
                    overLocs.add(loc)
        if len(overLocs)<1: overLocs=''
        return overLocs
    def nearGenes(self,loc,up=-5000,dn=0):
        tloc=Loc(loc)
        overLocs=set()
        if not tloc.chr in self.chrs:
            return overLocs
        for loc in self.chr2locs[tloc.chr]:
            aa,bb=loc.regionUpDn('TSS','TTS',up,dn)
            if tloc.la<=aa<=tloc.lb or tloc.la<=bb<=tloc.lb or aa<=tloc.la<=bb:
                overLocs.add(loc)
                continue
            #if tloc.la>loc.TSS: ###loc locate tloc's left
            #    lldis=tloc.la - loc.TSS
            #    if loc.sign=='+':
            #        ldis=dn
            #    else:
            #        ldis=up
            #    if lldis<=ldis:
            #        overLocs.add(loc)
            #else:              ### loc locate tloc's right
            #    rrdis=loc.TSS - tloc.lb
            #    #if rrdis<0:  ##check error
            #    #    print tloc,rloc,rloc.loc,'\n'
            #    if loc.sign=='+':
            #        rdis=up
            #    else:
            #        rdis=dn
            #    if rrdis<=rdis:
            #        overLocs.add(loc)
        if len(overLocs)<1: overLocs=''
        return overLocs

class RefSeq:
    def __init__(self,xl):
        self.geneid=''
        self.des=''
        self.sym,self.ref,self.chr,self.sign=xl[0:4]
        self.data=xl
        self.la,self.lb,self.lc,self.ld=map(int,xl[4:8])
        try:
            self.starts=map(int,xl[9].split(',')[:-1])
        except:
            print(xl)
        self.ends=map(int,xl[10].split(',')[:-1])
        self.geneLength=self.lb-self.la
        self.TSS,self.TTS=self.la,self.lb
        self.trS,self.trE=xl[6],xl[7]
        self.loc='%s:%s-%s:%s'%(self.chr,self.la,self.lb,self.sign)
        #self.loc='%s:%s:%s-%s:%s'%(self.ref,self.chr,self.la,self.lb,self.sign)
        if self.sign=='-': 
            self.TSS,self.TTS=self.lb,self.la
            self.trS,self.trE=self.trE,self.trS
        
    def getRefSeq2eles(self):    ### updated by sszhu1007@gmail.com 11.11.17 12:29:28 ###
        self.eles=Ref2eles(self.data)
    def getPairEles(self,tag1,tag2,locs):
        aa=int(self.__dict__[tag1])
        bb=int(self.__dict__[tag2])
        olocs=[]
        for a,b in locs:
            if self.sign=='+':
                olocs.append((aa+int(a),bb+int(b)))
            else:
                olocs.append((bb-int(b),aa-int(a)))
        return olocs
    def getEles(self,tag,locs):
        aa=int(self.__dict__[tag])
        olocs=[]
        for a,b in locs:
            if self.sign=='+':
                olocs.append((aa+int(a),aa+int(b)))
            else:
                olocs.append((aa-int(b),aa-int(a)))
        return olocs
    def getRegionEles(self,tag1,tag2,locs):
        #### return pair tuples as defined TSS,TTS,((-5000,0),)
        aa=int(self.__dict__[tag1])
        bb=int(self.__dict__[tag2])
        olocs=[]
        for a,b in locs:
            if self.sign=='+':
                olocs.append((aa+int(a),bb+int(b)))
            else:
                olocs.append((bb-int(b),aa-int(a)))
        return olocs
    def tagUpDn(self,tag,up,dn):
        ###up,dn=ref.tagUpDn('TSS',-10000,5000)
        aa=int(self.__dict__[tag])
        olocs=''
        if self.sign=='+':
            olocs=(aa+int(up),aa+int(dn))
        else:
            olocs=(aa-int(dn),aa-int(up))
        return olocs
    def regionUpDn(self,tag1,tag2,up,dn):
        aa=int(self.__dict__[tag1])
        bb=int(self.__dict__[tag2])
        olocs=''
        if self.sign=='+':
            olocs=(aa+int(up),bb+int(dn))
        else:
            olocs=(bb-int(dn),aa-int(up))
        return olocs

    ##def peakLoction(self,loc,up=1000,dn=1000):
    ####
        '''
        get loction for your peak, in check region
        '''
    def peakLocation(self,loc):
        tloc=Loc(loc)
        aa,bb=tloc.la,tloc.lb
        la,lb=self.la,self.lb
        if self.chr != tloc.chr:
            return 'out','diffChr'
        elif la<=aa<bb<=lb:
            return 'in',min(aa-la,lb-bb)
        elif aa<=la<lb<=bb:
            return 'contain',min(la-aa,bb-lb)
        elif la<=aa<=lb or la<=bb<=lb:
            return 'cross',0
        else:
            return 'out',min(abs(la-bb),abs(aa-lb))
    ####
    #def __str__(self):
    def __eq__(self,tt):    ### updated by sszhu1007@gmail.com 12.09.25 11:36:53 ###
        return self.loc == tt.loc
    def __hash__(self):
        return str.__hash__(self.loc)
    def __repr__(self):
        return self.sym
class Loc:
    def __init__(self,loc):
        self.loc=loc
        xl=re.split(':|\|',loc)
        self.chr=xl[0]
        #self.la,self.lb=re.search('([\d.]+)-([\d.]+)',xl[1]).groups()
        if re.search('^\d+$',xl[1]):
            self.la=int(xl[1])-1    ### updated by sszhu1007@gmail.com 12.04.11 17:49:24 ###
            self.lb=int(xl[1])
        elif re.search('(\d+)-(\d+)',xl[1]):
            self.la,self.lb=re.search('(\d+)-(\d+)',xl[1]).groups()
            self.la,self.lb=int(self.la),int(self.lb)
        if len(xl)>2: 
            self.sign=xl[2]
        else: self.sign='+'
    def __eq__(self):
        return self.loc
    def __str__(self):
        return self.loc
    def __hash__(self):
        return str.__hash__(self.loc)
class Bed:
    def __init__(self,file=os.path.join(dataDir,'test.bed')):
        self.chrs=set()
        self.locs=set()
        self.chr2locs={}
        for line in open(file):
            xl=line.rstrip().split('\t')
#import re
#re.search('(:.*-)','chr23:234234-234234').groups()
#re.search('^(chr\S+:\d+)$','chr23:234234').groups()
#mm=re.search('^(chr\S+):(\d+)$','chr23:234234').groups()
#tt,aa=re.search('^(chr\S+):(\d+)$','chr23:234234').groups()
            if re.search('^chr\S+:\d+-\d+', xl[0]):
                #chr,la,lb=re.search('^(chr\S+):(\d+)-(\d+)',xl[0])
                loc=xl[0]
            elif re.search('^chr\S+:\d+$',xl[0]):
                chr,la=re.search('^(chr\S+):(\d+)$',xl[0])
                loc='%s:%s-%s'%(chr,int(la)-1,la)
            elif len(xl)<3: continue
            #### use loc as unique id
            else:
                if len(xl)<6: sign=xl[3]
                else: sign=xl[5]
                loc='%s:%s-%s:%s'%(xl[0],xl[1],xl[2],sign)
            if loc in self.locs: continue
            self.locs.add(loc)
            obj=Loc(loc)
            self.chrs.add(obj.chr)
            self.chr2locs.setdefault(obj.chr,set()).add(obj)
    def siteNearbyLocs(self,site,cut=100000):
        chr,mid=site.split(":")
        mid=int(mid)
        ##ldis=rdis=sys.maxint
        ldis=rdis=100000000
        lloc=rloc='NA'
        overLocs=set()
        if not chr in self.chrs:
            return ldis,lloc,overLocs,rdis,rloc
        for loc in self.chr2locs[chr]:
            if loc.la<=mid<=loc.lb:
                overLocs.add(loc)
                continue
            if mid>loc.lb: ###loc locate tloc's left
                lldis=mid - loc.lb
                if lldis<ldis:
                    ldis=lldis
                    lloc=loc
            else:              ### loc locate tloc's right
                rrdis=loc.la - mid
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if rrdis<rdis:
                    rdis=rrdis
                    rloc=loc
        return ldis,lloc,overLocs,rdis,rloc
    def locNearbyLocs(self,loc):    ### updated by sszhu1007@gmail.com 12.04.28 14:06:40 ###  1 based inclusion
        tloc=Loc(loc)
        ##ldis=rdis=sys.maxint
        ldis=rdis=100000000
        lloc=rloc='NA'
        overLocs=set()
        if not tloc.chr in self.chrs:
            return ldis,lloc,overLocs,rdis,rloc
        for loc in self.chr2locs[tloc.chr]:
            if loc.la<=tloc.la<=loc.lb or loc.la<=tloc.lb<=loc.lb \
                or tloc.la<=loc.la<=tloc.lb:
                overLocs.add(loc)
                continue
            if tloc.la>loc.lb: ###loc locate tloc's left
                lldis=tloc.la - loc.lb
                if lldis<ldis:
                    ldis=lldis
                    lloc=loc
            else:              ### loc locate tloc's right
                rrdis=loc.la - tloc.lb
                #if rrdis<0:  ##check error
                #    print tloc,rloc,rloc.loc,'\n'
                if rrdis<rdis:
                    rdis=rrdis
                    rloc=loc
        return ldis,lloc,overLocs,rdis,rloc

if __name__ == '__main__':
    #Refs()
    from refFlat import RefFlat
    refFlat=RefFlat("/picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/refFlat.txt")
    refFlat.getSym2ref()
    refFlat.printSym2merge()
    #refFlat=RefFlat()
    #print refFlat.nearbyLocs('chr1:20000-400000:-')
    #print refFlat.nearbyTSS('chr2:95374953-95375271')
    #print Centromere().disFromCentr('chr1:20000-400000:-')
    #from refFlat import ref2infor
    #
    #ref2infor=ref2infor('refFlat.txt')
