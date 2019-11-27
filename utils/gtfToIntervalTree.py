# -*- author: Zheng Luo <crescentluozheng@gmail.com> -*-

import os
import cgranges as cr

anno = cr.cgranges()

anno.add("chr1", 10, 20, 0)
anno.add("chr1", 15, 25, 1)
anno.add("chr1", 13, 15, 1)
anno.index()
for st, en, label in anno.overlap("chr1", 12, 16):
        print(st, en, label)
print anno.coverage("chr1", 12, 16)
