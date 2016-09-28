#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

dir_data = "/Users/KANG/geneoscopy_dev/data/20160903_project_combined_2283_a_b_c_d/"
file_expr = dir_data + "chipdata_geneset_x_valid_chips.expression_console.txt"
file_sample = dir_data + "sample_sheet_combined_a_b_c_d.txt"
label = "P"

expr = np.loadtxt(file_expr, dtype=str, delimiter='\t')
meta = np.loadtxt(file_sample, dtype=str, usecols=[1,2,3], skiprows=1, delimiter='\t')

batch = {'1':[],'2':[],'3':[],'4':[]}
for i in range(len(meta)):
	if meta[i,2].split('.')[1] == label:
		tmp = meta[i,0].split('-')
		if len(tmp) == 1:
			batch['1'].append(".".join([meta[i,1],label]))
		else:
			batch[tmp[1]].append(".".join([meta[i,1],label]))

cg_genes = {'NDRG4':'TC16000503.hg.1', 'BMP3':'TC04000449.hg.1', 'KRAS':'TC12001314.hg.1'}
gene_indx = {}
gene_expr = {}
for g in cg_genes.keys():
	gene_indx[g] = np.where(expr[:,0] == cg_genes[g])[0][0]
	gene_expr[g] = []
	for k in sorted(batch.keys()):
		gene_expr[g].append([])
		for i in range(len(batch[k])):
			if batch[k][i] in expr[0,:]:
				sample_indx = np.where(expr[0,:] == batch[k][i])[0][0]
				gene_expr[g][len(gene_expr[g])-1].append(float(expr[gene_indx[g], sample_indx]))

for g in cg_genes.keys():
	plt.figure()
	bp = plt.boxplot(gene_expr[g])
	plt.title(g)
	plt.xticks(range(1,5),['a','b','c','d'])
	plt.xlabel('batch')
	plt.ylabel('log expression')
	plt.savefig(dir_data+"figures/"+g+"_"+label+".png", format="png")
