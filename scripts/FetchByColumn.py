"""Fetch by Colunmn"""
import sys
INDEX_SET = set()
with open(sys.argv[1]) as index_file:
    for line in index_file:
        sline = line.rstrip().split()
        INDEX_SET.add(sline[0])
with open(sys.argv[2]) as data_file:
    FILE_HEADER = data_file.readline()
    print FILE_HEADER.rstrip()
    for line in data_file:
        sline = line.rstrip().split()
        if sline[0] in INDEX_SET:
            print line.rstrip()
