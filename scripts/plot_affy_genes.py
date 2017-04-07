#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

file_expr = '/Users/KANG/geneoscopy_dev/data/20170217_combined_projects_batch1-17/chipdata_rma.expression_console.sst_rma.txt'
file_conv = '/Users/KANG/geneoscopy_dev/data/external_data/Genecards_colon_cancer/GeneCards_Nanostring_CIViC_genes_annotated.txt'
gene_list = ['PTGS2']

##parse expr data
affy_data = np.loadtxt(file_expr, dtype=str, delimiter="\t")
affy_genes = affy_data[1:,0]
affy_labels = affy_data[0,1:]
affy_data = np.array(affy_data[1:,1:], dtype=float)

indx_c = []
indx_p = []
indx_n = []
for i in range(len(affy_labels)):
	tmp =  affy_labels[i].split()[1].split('.')[0]
	if tmp == 'N':
		indx_n.append(i)
	elif tmp == 'P':
		indx_p.append(i)
	elif tmp == 'C':
		indx_c.append(i)

##parse conversion data
conv_dict = {}
f = open(file_conv, "r")
lines = f.readlines()
for i in range(len(lines)):
	line = lines[i].strip().split("\t")
	if len(line) >2:
		conv_dict[line[1]] = line[2].split(",")
f.close()


##make boxplot plot for each gene
colors= ["r", "b", "g"]
for g in gene_list:
	tcs = conv_dict[g]
	for tc in tcs:
		indx_g = np.where(affy_genes == tc)[0][0]
		tmp_data = [affy_data[indx_g,indx_c],
					affy_data[indx_g,indx_p],
					affy_data[indx_g,indx_n]]
		##boxplot
		plt.figure(num=None, figsize=(4,3), dpi=150)
		bp = plt.boxplot(tmp_data, 0, "", patch_artist=True) 
		for b,c in zip(bp['boxes'], colors):
			plt.setp(b, facecolor=c, alpha=.5)
		##overlay scatter points
		for i in range(3):
			plt.scatter(np.random.normal(i+1, .025, size=len(tmp_data[i])), tmp_data[i], 
						facecolors='none', linewidth=1.5)
		plt.xticks([1,2,3],["Cancer","Polyp","Normal"])
		plt.ylim([4.7,6.1])
		plt.title(g)
		plt.show()
