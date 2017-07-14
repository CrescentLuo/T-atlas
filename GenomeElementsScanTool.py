__author__ = 'luozheng'
import sys
import argparse
import intervaltree

geneSet = []
chromDict = {}
refFile = sys.argv[1]
regionFile = sys.argv[2]
print refFile
def readRefFile(refFile):
    with open(refFile) as ref:
        for line in ref:
            sline = line.rstrip().split()
            gene = sline[0]
            chrom = sline[2]
            strand = sline[3]
            txS = sline[4]
            txE = sline[5]
            cdsS = sline[6]
            cdsE = sline[7]
            exonCnt = int(sline[8])
            exonS = sline[9].rstrip(',').split(',')
            exonE = sline[10].rstrip(',').split(',')
            transcript = [gene, chrom, strand, txS, txE, cdsS, cdsE, exonCnt, exonS, exonE]
            geneSet.append(transcript)
            if chrom not in chromDict:
                chromDict[chrom] = len(chromDict)
                Exon.append([])
                Intron.append([])
                CDS.append([])
                UTR_5.append([])
                UTR_3.append([])
    geneSet.sort()
CDS = []
Exon = []
UTR_5 = []
UTR_3 = []
ExonIntersection = []
Intron = []
#rRNA = []
def CommonSegment(a,b,c,d):
    commonLen = 0
    if c > b:
        commonLen =  0
    elif c < b and c >= a:
        if d < b:
            commonLen =  d - c
        else:
            commonLen = b - c
    elif c < a:
        if d >= b:
            commonLen =  b - a
        elif d >= a and d < b:
            commonLen =  d - a
        else:
            commonLen =  0
    ALen = b - a +0.0
    BLen = d - c +0.0
    if commonLen /ALen >=0.05 or commonLen/BLen >= 0.05:
        return True
    return False
def getElements():
    Flag = True
    for gene in geneSet:
        if Flag:
            print gene
            Flag = False
        cdsS = int(gene[5])
        cdsE = int(gene[6])
        exonCnt = gene[7]
        exonS = gene[8]
        exonE = gene[9]
        chrom = gene[1]
        strand = gene[2]
        for i in range(exonCnt):
            exonS[i] = int(exonS[i])
            exonE[i] = int(exonE[i])

        if cdsS != cdsE:
            for i in range(exonCnt):
                if exonE[i] <= cdsS or exonS[i] >= cdsE:
                    if exonE[i] <= cdsS:
                        if strand == "+":
                            UTR_5[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
                        else:
                            UTR_3[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
                    elif exonS[i] >= cdsE:
                        if strand == "+":
                            UTR_3[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
                        else:
                            UTR_5[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
                    continue
                else:
                    if exonE[i] <= cdsE:
                        if exonS[i] >= cdsS:
                            CDS[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
                        else:
                            CDS[chromDict[chrom]].append([chrom,cdsS,exonE[i],gene[2],gene[0]])
                            if strand == "+":
                                UTR_5[chromDict[chrom]].append([chrom,exonS[i],cdsS,gene[2],gene[0]])
                            else:
                                UTR_3[chromDict[chrom]].append([chrom,exonS[i],cdsS,gene[2],gene[0]])
                    else:
                        if exonS[i] >= cdsS:
                            CDS[chromDict[chrom]].append([chrom,exonS[i],cdsE,gene[2],gene[0]])
                            if strand == "+":
                                UTR_3[chromDict[chrom]].append([chrom,cdsE,exonE[i],gene[2],gene[0]])
                            else:
                                UTR_5[chromDict[chrom]].append([chrom,cdsE,exonE[i],gene[2],gene[0]])
                        else:
                            CDS[chromDict[chrom]].append([chrom,cdsS,cdsE,gene[2],gene[0]])
                            if strand == "+":
                                UTR_5[chromDict[chrom]].append([chrom,exonS[i],cdsS,gene[2],gene[0]])
                                UTR_3[chromDict[chrom]].append([chrom,cdsE,exonE[i],gene[2],gene[0]])
                            else:
                                UTR_3[chromDict[chrom]].append([chrom,exonS[i],cdsS,gene[2],gene[0]])
                                UTR_5[chromDict[chrom]].append([chrom,cdsE,exonE[i],gene[2],gene[0]])
        for i in range(exonCnt):
            Exon[chromDict[chrom]].append([chrom,exonS[i],exonE[i],gene[2],gene[0]])
            if i:
                Intron[chromDict[chrom]].append([chrom,exonE[i-1],exonS[i],gene[2],gene[0]])
def ElementSort():
    for chrom in Exon:
        chrom.sort()
    for chrom in Intron:
        chrom.sort()
    for chrom in CDS:
        chrom.sort()
    for chrom in UTR_5:
        chrom.sort()
    for chrom in UTR_3:
        chrom.sort()
readRefFile(refFile)
getElements()
ElementSort()
def biSearch(target,searchSet):
    left = 0
    right = len(searchSet) - 1
    while(left < right):
        mid = (left + right) /2
        #print left,right,mid
        midRegion = searchSet[mid]
        if not(target[0] >= midRegion[2] or target[1] <= midRegion[1]) and CommonSegment(target[0],target[1],midRegion[1],midRegion[2]):
                return mid
        else:
            if target[0] >= midRegion[1]:
                left = mid + 1
            else:
                right = mid
    if left == right:
        midRegion = searchSet[left]
        if CommonSegment(target[0],target[1],midRegion[1],midRegion[2]):
            return left
    return -1
InterGenicRegionCnt = 0
ExonRegionCnt = 0
IntronRegionCnt = 0
CDSRegionCnt = 0
UTR_5_RegionCnt = 0
UTR_3_RegionCnt = 0
CodingReigonCnt = 0
NonCodingRegionCnt = 0
AlternativedRegionCnt = 0
rRNARegionCnt = 0
snoRNARegionCnt = 0
#rRNA[chromDict["chrUn_gl000220"]].append(["chrUn_gl000220",105424,105424,"+","rRNA"])
with open(regionFile) as rf:
    rf.readline()
    for line in rf:
        sline = line.rstrip().split()
        chrom = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        if chrom in chromDict:
            chrN = chromDict[chrom]
        else:
            continue
        exonNum = len(Exon[chromDict[chrom]])
        ExonFlag = False
        IntronFlag = False
        CDSFlag = False
        UTR_5Flag = False
        UTR_3Flag = False
        intronNum = len(Intron[chromDict[chrom]])
        if (biSearch([start,end],Exon[chrN])) != -1:
            ExonFlag = True
            if Exon[chrN][biSearch([start,end],Exon[chrN])][4]=="RNA5-8S5":
                rRNARegionCnt +=1
            geneName = Exon[chrN][biSearch([start,end],Exon[chrN])][4]
            if len(geneName) >= 6  and (geneName[0:5]=="SNORD" or geneName[0:5]=="SNORA" or geneName[0:6]=="SCARNA"):
                snoRNARegionCnt += 1
            ExonRegionCnt +=1
        if (biSearch([start,end],Intron[chrN])) != -1:
            IntronFlag =True
            IntronRegionCnt +=1
        if (biSearch([start,end],CDS[chrN])) != -1:
            CDSFlag = True
            CDSRegionCnt +=1
        if (biSearch([start,end],UTR_5[chrN])) != -1:
            UTR_5Flag = True
            UTR_5_RegionCnt +=1
        if (biSearch([start,end],UTR_3[chrN])) != -1:
            UTR_3Flag = True
            UTR_3_RegionCnt +=1
        if ExonFlag and IntronFlag:
            AlternativedRegionCnt += 1
        if not(ExonFlag) and not(IntronFlag):
            InterGenicRegionCnt += 1
        if (UTR_3Flag or UTR_5Flag or CDSFlag):
            CodingReigonCnt +=1
        elif ExonFlag:
            NonCodingRegionCnt +=1

    print "ExoRegionCnt:",ExonRegionCnt
    print "IntronRegionCnt:",IntronRegionCnt
    print "InterGenicRegionCnt",InterGenicRegionCnt
    print "AlternativedRegionCnt",AlternativedRegionCnt
    print "CDSRegionCnt",CDSRegionCnt
    print "5UTRRegionCnt",UTR_5_RegionCnt
    print "3UTRRegionCnt",UTR_3_RegionCnt
    print "CodingRegionCnt",CodingReigonCnt
    print "NonCodingRegionCnt",NonCodingRegionCnt
    print "rRNARegionCnt",rRNARegionCnt
    print "snoRNARegionCnt",snoRNARegionCnt
