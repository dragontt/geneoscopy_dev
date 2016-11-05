#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

dir_data = "/Users/KANG/geneoscopy_dev/data/20161014_project_combined_2283_abcde/"
file_expr = dir_data + "batch_effect_removal/chipdata_geneset_x_valid_chips.batcheff_removed.txt"
file_sample = dir_data + "sample_sheet_combined_abcde.two_group_no_benign.txt"

bs = ["1", "2", "3", "4", "5"]
labels = ["C", "N"]
colors_arr = ["r", "g"]
"""
# cg_genes = {'NDRG4':'TC16000503.hg.1', 'BMP3':'TC04000449.hg.1', 'KRAS':'TC12001314.hg.1'}
cg_genes = {'TC6_ssto_hap7000119' : 'TC6_ssto_hap7000119.hg.1'}
"""
cg_genes = {}
file_genes = dir_data + "../external_data/CIViC/civic_genes_TCs.txt"
lines = open(file_genes, "r").readlines()
for line in lines:
	tmp = line.strip().split('\t')
	if len(tmp) > 1:
		gene = tmp[0]
		tcs = tmp[1].split(",")
		for tc in tcs:
			cg_genes[gene +"_"+ tc] = tc

expr = np.loadtxt(file_expr, dtype=str, delimiter='\t')
meta = np.loadtxt(file_sample, dtype=str, usecols=[1,2,3], skiprows=1, delimiter='\t')

## identify batches from sample sheet
batch = {}
for i in range(len(bs)):
	for j in range(len(labels)):
		batch[bs[i]+labels[j]] = []

for i in range(len(labels)):
	label = labels[i]
	for i in range(len(meta)):
		if meta[i,2].split('.')[1] == label:
			tmp = meta[i,0].split('-')
			if len(tmp) == 1:
				batch['1'+label].append(".".join([meta[i,1],label]))
			else:
				batch[tmp[1]+label].append(".".join([meta[i,1],label]))

## get expression for each batch
gene_indx = {}
gene_expr = {}
for g in cg_genes.keys():
	gene_indx[g] = np.where(expr[:,0] == cg_genes[g])[0][0]
	gene_indx[g] = gene_indx[g]
	gene_expr[g] = []
	for k in sorted(batch.keys()):
		gene_expr[g].append([])
		for i in range(len(batch[k])):
			if batch[k][i] in expr[0,:]:
				sample_indx = np.where(expr[0,:] == batch[k][i])[0][0]
				gene_expr[g][len(gene_expr[g])-1].append(float(expr[gene_indx[g], sample_indx]))
	## add all samples expressions for each batch
	all_expr = []
	for j in range(len(labels)):
		all_expr.append([])
	for i in range(len(bs)):
		for j in range(len(labels)):
			all_expr[j] += gene_expr[g][i*len(labels)+j]
	gene_expr[g] += all_expr

## plot boxplot
bs.append("all")
locs = []
for i in range(len(bs)):
	locs += [x+i*(len(labels)+1) for x in range(len(labels))]
colors = colors_arr*(len(bs))
for i in range(len(labels)):
	batch["all"+labels[i]] = []

for g in cg_genes.keys():
	plt.figure()
	bp = plt.boxplot(gene_expr[g], positions=locs, patch_artist=True)
	for b,c in zip(bp['boxes'], colors):
		plt.setp(b, facecolor=c)
	plt.title(g)
	plt.xticks(locs, sorted(batch.keys()))
	plt.xlabel('batch + label')
	plt.ylabel('log expression')
	plt.savefig(dir_data+"figures/"+g+".png", format="png")

	# print g, "\t", len(bp["fliers"][21].get_data()[1]) ,"\t", len(bp["fliers"][23].get_data()[1])


