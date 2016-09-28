#!/usr/bin/python
import sys
import numpy as np

def ParseGenes(transcripts, genes, delim=None):
	if delim is None:
		delim = ","
	unique_transcripts = []
	gene_dict = {}
	for i in range(len(transcripts)):
		if transcripts[i] != '':
			ts = transcripts[i].split(delim)
			for t in ts:
				if t != '':
					t = t.strip('"').strip().split(".")[0]
					unique_transcripts.append(t)
					gene_dict[t] = genes[i]
	return (np.array(unique_transcripts), gene_dict)

def CompareGenes(g_genes, e_genes, e_data, gene_dict):
	intersected_genes = []
	# print np.intersect1d(g_genes, e_genes)
	for x in np.intersect1d(g_genes, e_genes):
		if gene_dict:
			intersected_genes.append((x, gene_dict[x]))
		else:
			intersected_genes.append(x)
	return intersected_genes

def LoadTable(filename, usecols, skiprows):
	f = open(filename, "r")
	lines = f.readlines()[0].split("\r")
	data = []
	for i in range(skiprows, len(lines)):
		line = lines[i].split("\t")
		tmp = []
		print i, line[1], line[4], line[6]
		for j in range(len(usecols)):
			tmp.append(line[usecols[j]])
		data.append(tmp)
	f.close()
	return np.array(data)

# read data
filename_geneoscopy = sys.argv[1]
filename_external_data = sys.argv[2]
# geneoscopy = np.loadtxt(filename_geneoscopy, dtype=str, delimiter='\t', usecols=[1,4,6], skiprows=1)
# delim = ","
geneoscopy = np.loadtxt(filename_geneoscopy, dtype=str, delimiter='\t', usecols=[2,1,3], skiprows=1)
delim = ","

external_data = np.loadtxt(filename_external_data, dtype=str, delimiter='\t')

# compare transcripts 
(geneoscopy_transcripts, foo) = ParseGenes(geneoscopy[:,0], geneoscopy[:,2], delim)
(external_data_transcripts, external_data_gene_dict) = ParseGenes(external_data[:,2], external_data[:,0])
intersected_transcripts = CompareGenes(geneoscopy_transcripts, external_data_transcripts, external_data, external_data_gene_dict)
print "Supported ENST transcripts", intersected_transcripts

# compare gene symbol
geneoscopy_symbol = geneoscopy[:,1]
for i in range(len(geneoscopy_symbol)):
	geneoscopy_symbol[i] = geneoscopy_symbol[i].strip()
external_data_symbol = external_data[:,0]
intersected_symbol = CompareGenes(geneoscopy_symbol, external_data_symbol, external_data, None)
print "Supported gene symbols", intersected_symbol

# compare accession
if external_data.shape[1] > 3:
	(geneoscopy_accession, foo) = ParseGenes(geneoscopy[:,2], geneoscopy[:,2], delim)
	(external_data_accession, external_data_gene_dict) = ParseGenes(external_data[:,3], external_data[:,0])
	intersected_symbol = CompareGenes(geneoscopy_accession, external_data_accession, external_data, external_data_gene_dict)
	print "Supported mRNA acession", intersected_symbol

