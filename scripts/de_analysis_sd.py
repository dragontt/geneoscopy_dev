#!/usr/bin/python
import numpy as np
import sys

# dir_proj = "/Users/KANG/geneoscopy_dev/data/run_proj_a_b_c_d/"
# file_expr_data = dir_proj + "training/chipdata.txt"
# file_valid_chips = dir_proj + "training/valid_chips.txt"
# file_de_genes = dir_proj + "training/top_de_genes_sd.txt"

file_expr_data = sys.argv[1]
file_valid_chips = sys.argv[2]
file_de_genes = sys.argv[3]

expr_data = np.loadtxt(file_expr_data, dtype=str, delimiter="\t")
genes = expr_data[1:,0]
samples = expr_data[0,1:]
expr_data = np.exp2(np.array(expr_data[1:,1:], dtype=float))
valid_chips = np.loadtxt(file_valid_chips, dtype=str)

label_dict = {}
for i in range(len(valid_chips)):
	label_dict[valid_chips[i,0]] = valid_chips[i,1]

indx_n = []
indx_c = []
for i in range(len(samples)):
	if label_dict[samples[i]] == '1':
		indx_c.append(i)
	else:
		indx_n.append(i)
print "Normal samples : cancer samples =", len(indx_n), ":", len(indx_c)

de_genes = []
for i in range(len(genes)):
	mean_n = np.mean(expr_data[i,indx_n])
	mean_c = np.mean(expr_data[i,indx_c])
	std_n = np.std(expr_data[i,indx_n])
	std_c = np.std(expr_data[i,indx_c])
	if (std_c > std_n) and (abs(std_c-std_n) > .1):
		std_c_counts = 0
		std_c_counts += len(np.where(expr_data[i,indx_c] > mean_c+1.5*std_c)[0])
		std_c_counts += len(np.where(expr_data[i,indx_c] < mean_c-1.5*std_c)[0])
		std_c_pct = std_c_counts/float(len(indx_c))
		if (std_c_pct > .1) and (abs(std_c-std_n) > .8):
			de_genes.append([genes[i], mean_n, std_n, mean_c, std_c, std_c_pct])

print "DE genes identified:", len(de_genes)
de_genes_header = np.array(['TC id', 'normal mean', 'normal std', 'cancer mean', 'cancer std', 'pct cancer > 1.5std'])
de_genes = np.vstack((de_genes_header[np.newaxis], np.array(de_genes, dtype=str)))
np.savetxt(file_de_genes, de_genes, fmt="%s", delimiter="\t")

