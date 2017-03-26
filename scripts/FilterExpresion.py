import sys
group_dict = {}
group_cnt = {}
exp_thres = float(sys.argv[3])
with open(sys.argv[1]) as group_info:
    header = group_info.readline()
    for line in group_info:
        sline = line.split()
        group_dict[sline[0]] = sline[1]
        if sline[1] not in group_cnt:
            group_cnt[sline[1]] = 1
        else:
            group_cnt[sline[1]] += 1 
with open(sys.argv[2]) as exp_matrix, open("Group_FPKM_ave.txt", "w") as group_out, open("Filtered_Res.txt", "w") as res_file:
    header = exp_matrix.readline()
    sample_info = header.rstrip().split()[1:-1]
    res_file.write(header)
    gene_ave = {}
    for sample in group_dict:
        gene_ave[group_dict[sample]] = 0.0
    group = [str(key) for key in group_cnt]
    header_line = "gene\t"+'\t'.join(group)+'\n'
    group_out.write(header_line)
    for line in exp_matrix:
        sline = line.rstrip().split()
        gsymbol = sline[0]
        gexp = [float(exp) for exp in sline[1:-1]]
        for i in range(len(sample_info)):
            group = group_dict[sample_info[i]]
            gene_ave[group] += gexp[i]
        flag = False
        for group in gene_ave:
            gene_ave[group] = gene_ave[group] / group_cnt[group]
            if gene_ave[group] >= exp_thres:
                flag = True
        if flag:
            ave_values = [gene_ave[group] for group in gene_ave]
            ave_values = [str(i) for i in ave_values]
            group_out.write(gsymbol+'\t'+'\t'.join(ave_values)+'\n')
            res_file.write(line)
