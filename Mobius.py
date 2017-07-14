#%%

import sys
import numpy as np 
import pandas as pd 
import seaborn as sns
from IPython.display import Image

class circPos:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

class circRNA:
    """circular RNA from CIRCexplorer2"""
    def __init__(self, circInfo):
        #print circInfo
        self.circPos = circPos(circInfo[0],int(circInfo[1]), int(circInfo[2]))
        self.name = circInfo[3]
        self.score = circInfo[4]
        self.strand = circInfo[5]
        self.thickStart = circInfo[6] 
        self.thickEnd = circInfo[7]
        self.itemRgb = circInfo[8]
        self.exonCount = circInfo[9]
        self.exonSizes = circInfo [10]
        self.exonOffsets = circInfo[11]
        self.readNumber = circInfo[12]
        self.circType = circInfo[13]
        self.geneName = circInfo[14]
        self.isoformName = circInfo[15]
        self.index =  circInfo[16]
        self.flankIntron = circInfo[17]
    def getReadCount(self) :
        return self.readNumber

class circSet:
    """circular RNA list"""
    def __init__(self,totalMapped):
        self.circSet = []
        self.totalMapped = totalMapped
    def addCirc(self, circRNA):
        self.circSet.append(circRNA)
    def getCircRNACount(self):
        return len(self.circSet)
    def getDist(self):
        readNum = {}
        for circ in self.circSet:
            if circ.readNumber in readNum:
                readNum[circ.readNumber] += 1
            else:
                readNum[circ.readNumber] = 1
        readNum = pd.DataFrame(readNum.items(), columns=['readCount', 'circNum'])
        sns.set_style("ticks")
        Image(sns.barplot(readNum['readCount'], readNum['circNum']))
        sns.despine()
    def getRPM(self):
        circExpHigh = 0
        circExpLow = 0
        for circ in self.circSet:
            circRPM = float(circ.readNumber * 10**6) / self.totalMapped
            if circRPM >= 0.1:
                circExpHigh += 1
            else:
                circExpLow += 1
        print circExpHigh,circExpLow


def readCircTable(annotate_file,totalMapped):
    circBed = pd.read_table(annotate_file, header=0)
    print circBed.shape
    circS = circSet(totalMapped)
    for idx in range(circBed.shape[0]):
        circ = circRNA(circBed.iloc[idx,:])
        circS.addCirc(circ)
    circS.getDist()
    return circS

#def plotCircCounts(*circSet):
    
    

circ_N_U_1 = readCircTable("/Users/luozheng/rnomics/Cooperation/Jeremy/DoG/1_CircRNAs/N_U_1.circ.txt",208433186)
circ_N_U_2 = readCircTable("/Users/luozheng/rnomics/Cooperation/Jeremy/DoG/1_CircRNAs/N_U_2.circ.txt",148979665)
circ_N_K_1 = readCircTable("/Users/luozheng/rnomics/Cooperation/Jeremy/DoG/1_CircRNAs/N_K_1.circ.txt",175877717)
circ_N_K_2 = readCircTable("/Users/luozheng/rnomics/Cooperation/Jeremy/DoG/1_CircRNAs/N_K_2.circ.txt",173577550)

cnt_wt_1 = circ_N_U_1.getRPM()
cnt_wt_2 = circ_N_U_2.getRPM()
cnt_kcl_1 = circ_N_K_1.getRPM()
cnt_kcl_2 = circ_N_K_2.getRPM()

count = pd.DataFrame([(cnt_wt_1,"wt","rep1"),(cnt_wt_2,"wt","rep2"),(cnt_kcl_1,"kcl","rep1"),(cnt_kcl_2,"kcl","rep2")], columns=['cnt', 'treatment','group'])
sns.barplot(x="treatment",y="cnt",hue="group",data=count)
