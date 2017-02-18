#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

file_expr = '/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_lit.txt'
gene_list = ['AKT1', 'AKT2', 'AKT3', 'MAK', 'BRCA1', 'BRCA2', 'ERBB2', 'ERBB4',
				'VEGF', 'VEGA', 'IL8', 'NFKB', 'IL22', 'PDL1', 'PDL', 'CTLA4']
dir_figures = '/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_figures/'

##parse data
ns_data = np.loadtxt(file_expr, dtype=str, delimiter="\t")
ns_genes = ns_data[2:,1]
ns_labels = ns_data[1,6:]
ns_data = np.array(ns_data[2:,6:], dtype=float)

indx_c = np.where(ns_labels == "Cancer")[0]
indx_p = np.where(ns_labels == "Polyp")[0]
indx_n = np.where(ns_labels == "Normal")[0]

##find avail and unavail Nanostring genes
genes_avail = []
genes_unavail = []
for g in gene_list:
	found = False
	for nsg in ns_genes:
		if nsg.startswith(g):
			genes_avail.append(nsg)
			found =True
	if not found:
		genes_unavail.append(g)
print "Available genes:", genes_avail
print "Unavailable genes:", genes_avail

##make boxplot plot for each gene
colors= ["r", "g", "b"]
for g in genes_avail:
	indx_g = np.where(ns_genes == g)[0][0]
	plt.figure(num=None, figsize=(5,4), dpi=150)
	tmp_data = [ns_data[indx_g,indx_c], 
				ns_data[indx_g,indx_p], 
				ns_data[indx_g,indx_n]]
	bp = plt.boxplot(tmp_data, 0, '', patch_artist=True) 
	for b,c in zip(bp['boxes'], colors):
		plt.setp(b, facecolor=c, alpha=.65)

	##overlay scatter points
	for i in range(3):
		plt.scatter(np.random.normal(i+1, .025, size=len(tmp_data[i])), tmp_data[i], 
					facecolors='none', linewidth=2)
	plt.xticks([1,2,3],["C","P","N"])
	plt.title(g)
	# plt.show()
	# plt.savefig(dir_figures + g + ".pdf", format="pdf")
