#%%
import sys
import argparse
import pandas as pd
scr = "/Users/luozheng/rnomics/Cooperation/Jeremy/CPSF73KD/PA1_Scr_circ_fusion.xlsx"
kd = "/Users/luozheng/rnomics/Cooperation/Jeremy/CPSF73KD/PA1_CPSF73KD_circ_fusion.xlsx"
scr = pd.read_excel(scr)
kd = pd.read_excel(kd)
scr['tag'] = scr['chrom'].map(str) +'@'+ scr['start'].map(str) +'@'+ scr['end'].map(str)
kd['tag'] = kd['chrom'].map(str) +'@'+ kd['start'].map(str) +'@'+ kd['end'].map(str)

kd = kd.merge(scr,how='left',on='tag')

kd.to_csv("/Users/luozheng/rnomics/Cooperation/Jeremy/CPSF73KD/merge.csv",index=False)  