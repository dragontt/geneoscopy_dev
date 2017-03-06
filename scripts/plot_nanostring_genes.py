#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

file_expr = '/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_lit.txt'
# gene_list = ['AKT1', 'AKT2', 'AKT3', 'MAK', 'BRCA1', 'BRCA2', 'ERBB2', 'ERBB4',
# 				'VEGF', 'VEGA', 'IL8', 'NFKB', 'IL22', 'PDL1', 'PDL', 'CTLA4']
gene_list = ['AKT2', 'IL22RA', 'NFKB1', 'VEGFA']
# gene_list = ['AKT2']
dir_figures = '/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_figures/'

##parse data
ns_data = np.loadtxt(file_expr, dtype=str, delimiter="\t")
ns_genes = ns_data[2:,1]
ns_labels = ns_data[1,6:]
ns_data = np.array(ns_data[2:,6:], dtype=float)

indx_c = np.where(ns_labels == "Cancer")[0]
indx_p = np.where(ns_labels == "Polyp")[0]
indx_n = np.where(ns_labels == "Normal")[0]

baseline_data = np.mean(ns_data[:,indx_n], axis=1)

##find avail and unavail Nanostring genes
genes_avail = []
# genes_unavail = []
for g in gene_list:
	found = False
	for nsg in ns_genes:
		if nsg.startswith(g):
			genes_avail.append(nsg)
			found =True
	# if not found:
	# 	genes_unavail.append(g)
print "Available genes:", genes_avail
# print "Unavailable genes:", genes_avail

##make boxplot plot for each gene
colors= [[35/255., 183/255., 176/255.], [0/255., 64/255., 100/255.]]
for g in genes_avail:
	indx_g = np.where(ns_genes == g)[0][0]
	plt.figure(num=None, figsize=(4,3), dpi=150)
	tmp_data = [ns_data[indx_g,np.append(indx_c,indx_p)],
				ns_data[indx_g,indx_n]]
	for i in range(len(tmp_data)):
		tmp_data[i] = np.log2(np.divide(tmp_data[i], baseline_data[indx_g]))
	##boxplot
	bp = plt.boxplot(tmp_data, patch_artist=True) 
	for b,c in zip(bp['boxes'], colors):
		plt.setp(b, facecolor=c)
	##overlay scatter points
	# for i in range(2):
	# 	plt.scatter(np.random.normal(i+1, .025, size=len(tmp_data[i])), tmp_data[i], 
	# 				facecolors='none', linewidth=1.5)
	plt.xticks([1,2],["Cancer/Adenomas","Healthy"])
	plt.title(g)
	# plt.show()
	plt.savefig(dir_figures + g + ".pdf", format="pdf")
