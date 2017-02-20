import sys
index_set = set()
with open(sys.argv[1]) as index_file:
    for line in index_file:
        sline = line.rstrip().split()
        index_set.add(sline[0])
with open(sys.argv[2]) as data_file:
    file_header = data_file.readline()
    print file_header.rstrip()
    for line in data_file:
        sline = line.rstrip().split()
        if sline[0] in index_set:
            print line.rstrip()
