import sys
import numpy
if len(sys.argv) == 5:
    col_s = int(sys.argv[3]) - 1
    col_e = int(sys.argv[4])
if len(sys.argv) >= 3:
    Outfile = sys.argv[2]
else:
    Outfile = "Zscore.txt"
with open(sys.argv[1]) as InputFile, open(Outfile,"w") as OutputFile:
    header = InputFile.readline()
    if len(sys.argv) == 5:
        header = header.split()
        header = [header[0]] + header[col_s:col_e]
        #print header
        header = '\t'.join(header)
    OutputFile.write(header)
    for line in InputFile:
        sline = line.split()
        if len(sys.argv) == 5:
            values = numpy.array(sline[col_s:col_e])
        else:
            values = numpy.array(sline[1:])
        values = values.astype(float)
        std = numpy.std(values)
        avg = numpy.mean(values)
        outline = []
        outline.append(sline[0])
        if avg ==0 and std ==0:
            continue
        for i in values:
            zscore = (i-avg)/std
            outline.append(str(zscore))
        OutputFile.write('\t'.join(outline)+'\n')
