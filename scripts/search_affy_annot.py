#!/usr/bin/python
import numpy as np

proj_dir = "/Users/KANG/geneoscopy_dev/data/"

file_annot = proj_dir + "HTA-2_0.na36.hg19.transcript.csv/HTA-2_0.na36.hg19.transcript.csv"
annot = np.loadtxt(file_annot, dtype=str, usecols=[0,7], delimiter=",")
"""
file_gene_list = proj_dir + "external_data/nanostring/PanCancer_nanostring_genes.txt"
file_gene_out = proj_dir + "external_data/nanostring/PanCancer_nanostring_genes_annotated.txt"
genes = np.loadtxt(file_gene_list, dtype=str)
for i in range(len(genes)):
	genes[i,0] = genes[i,0].split(".")[0]
"""
file_gene_list = proj_dir + "external_data/Genecards_colon_cancer/GeneCards_genes.txt"
file_gene_out = proj_dir + "external_data/Genecards_colon_cancer/GeneCards_genes_annotated.txt"
genes = np.loadtxt(file_gene_list, dtype=str)[np.newaxis].T

tc = []
for i in range(len(genes)):
	tc.append(set())

for i in range(len(annot)):
	if i % 1000 == 0:
		print "analyzing row", i
	for d in annot[i,1].strip('"').split('//'):
		d = d.strip()
		for j in range(len(genes)):
			for g in genes[j]:
				if d == g:
					tc[j].add(annot[i,0].strip('"'))

for i in range(len(tc)):
	tc[i] = ",".join(tc[i]) if tc[i] else ""
tc = np.array(tc)

genes = np.hstack((genes, tc[np.newaxis].T))
np.savetxt(file_gene_out, genes, fmt="%s", delimiter="\t")

