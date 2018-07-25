import sys

class GTF:
    """ GTF annotation """

    def __init__(self):
        self.chrom = '' # Chromosome
        self.source = '' # Source
        self.feature= '' # Feature type e.g gene/exon
        self.start = 0 #Feature start position
        self.end = 0 #Feature end position
        self.score = 0 # Feature score
        self.strand = '+'# Feature strand
        self.frame = ''
        self.attributes = '' # tag-value pairs
        self.attr_dict = {}
        self.is_exon = False # Ture if this feature is an exon
        self.attr_dict = dict()

    def parse_attribute(self):
        """ parse attibutes into dictionary"""
        self.attr_dict = dict([(i.split()[0],i.split()[1].strip('"')) \
                            for i in self.GTF.attributes.rstrip(';').split(';')])
    def get_attribute(self, attr):
        if attr in self.attr_dict:
            return self.attr_dict[attr]
        else:
            return "NA"


class GTFParser:
    """ GTF Parser """

    def __init__(self):

    def __init__(self,file):

    def __init__(self,line):
        sline = line.rstrip().split()
        self.GTF = GTF()
        if sline[2] == 'exon':
            self.GTF.is_exon = True
        self.GTF.chrom = sline[0]
        self.GTF.source = sline[1]
        self.GTF.feature = sline[2]
        self.GTF.start = int(sline[3])
        self.GTF.end = int(sline[4])
        self.GTF.score = sline[5]
        self.GTF.strand = sline[6]
        self.GTF.frame = sline[7][0]
        self.GTF.attributes = sline[8]
        self.transcript_dict = dict()
        return GTF
        
    
    def get_attribute(self):

    def add_exon_to_transcript(self, ts, exon):
        vector<string> attributes;
    Tokenize(gtf1.attributes, attributes, ';');
    string transcript_id = parse_attribute(attributes, "transcript_id");
    string gene_name = parse_attribute(attributes, "gene_name");
    //create a BED6 object
    BED exon = BED(gtf1.seqname, gtf1.start,
                   gtf1.end, gtf1.feature,
                   gtf1.score, gtf1.strand);
    if(transcript_id != string("NA")) {
        transcript_map_[transcript_id].exons.push_back(exon);
        set_transcript_gene(transcript_id, gene_name);
    }
}


    

    def sort_exon_within_transcript(self):
    
    def annotate_transcript_with_bins(self):

    def print_transcript(self):

    def transcript_from_bin(self):
    
    def ts_get_exon(self):
    def ts_get_gene(self):
    def ts_set_gene(self):    
       