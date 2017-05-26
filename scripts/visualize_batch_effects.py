#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dir_data = "/Users/KANG/geneoscopy_dev/data/20170217_combined_projects_batch1-17/"
file_expr = dir_data + "chipdata_geneset_x_valid_chips_full.txt"
# file_expr = dir_data + "chipdata_geneset_x_valid_chips_full.batch_removed.txt"
file_batch = dir_data + "sample_batch.txt"
dir_figures = "/Users/KANG/Desktop/tmp/"

labels = ["C", "N", "P"]
colors_arr = ["r", "g", "b"]
# """
# cg_genes = {'PTGS2': 'TC01003638.hg.1'}

# cg_genes = {'COMP':'TC19001293.hg.1', 
# 			'TAC1':'TC07002473.hg.1',
# 			'LRP5':'TC11000707.hg.1',
# 			'PTGES2':'TC09001623.hg.1',
# 			'HOXA10':'TC07002854.hg.1',
# 			'THBD':'TC20000702.hg.1',
# 			'EPOR':'TC19002676.hg.1',
# 			'BVES':'TC06001976.hg.1',
# 			'MUC6':'TC11001256.hg.1'}

# cg_genes = {'COMP':'TC19001293.hg.1',
# 			'THBD':'TC20000702.hg.1',
# 			'TAC1':'TC07002473.hg.1',
# 			'EPOR':'TC19002676.hg.1',
# 			'ENPP2':'TC08001557.hg.1',
# 			'PTGES2':'TC09001623.hg.1',
# 			'LRP5':'TC11000707.hg.1',
# 			'OLFM2':'TC19001155.hg.1',
# 			'CSNK2A1':'TC20001369.hg.1',
# 			'BVES': 'TC06001976.hg.1'}

cg_genes = {'KRAS': 'TC12001314.hg.1',
			'APC': 'TC05003398.hg.1'}

# """

"""
cg_genes = {}
file_conv = dir_data + "../external_data/Genecards_colon_cancer/GeneCards_Nanostring_CIViC_genes_annotated.txt"
conv_dict = {}
f = open(file_conv, "r")
lines = f.readlines()
for i in range(len(lines)):
	line = lines[i].strip().split("\t")
	if len(line) >2:
		gene = line[1]
		for tc in line[2].split(","):
			cg_genes[gene +"_"+ tc] = tc
f.close()
"""

expr = np.loadtxt(file_expr, dtype=str, delimiter='\t')
header = expr[0,1:]
meta = np.loadtxt(file_batch, dtype=str, skiprows=1, delimiter='\t')
meta_dict = {}
for i in range(len(meta)):
	meta_dict[meta[i,0]] = meta[i,1]
bs = list(np.array(sorted(np.unique(np.array(meta[:,1], dtype=int))), dtype=str))
print bs

## identify batches from sample sheet
batch = {}
for i in range(len(bs)):
	for j in range(len(labels)):
		batch[bs[i]+labels[j]] = []

for i in range(len(labels)):
	label = labels[i]
	for i in range(len(header)):
		tmp_id, tmp_label = header[i].split(".")
		if tmp_label == label:
			batch[meta_dict[tmp_id]+label].append(header[i])

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
os.makedirs(dir_figures)
bs.append("all")
locs = []
for i in range(len(bs)):
	locs += [x+i*(len(labels)+1) for x in range(len(labels))]
for i in range(len(labels)):
	batch["all"+labels[i]] = []

for g in cg_genes.keys():
	print "Plotting", g
	plt.figure(num=None, figsize=(20, 6), dpi=80)
	bp = plt.boxplot(gene_expr[g], positions=locs, patch_artist=True) 
	
	colors = colors_arr*(len(bs))
	empty_indx = [i for i in range(len(gene_expr[g])) if gene_expr[g][i] == []]
	colors = np.delete(colors, empty_indx)
	for b,c in zip(bp['boxes'], colors):
		plt.setp(b, facecolor=c)
	
	plt.title(g)
	plt.xticks(locs, sorted(batch.keys()))
	plt.xlabel('batch + label')
	plt.ylabel('log expression')
	# plt.show()
	plt.savefig(dir_figures + g + ".png", format="png")
