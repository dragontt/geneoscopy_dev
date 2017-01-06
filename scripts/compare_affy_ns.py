#!/usr/bin/python
import sys
import argparse
import numpy as np
# from scipy.stats.mstats import gmean
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Compare how Affy microarray data correlates to Nanostring nCounter data.')
	parser.add_argument('-m', '--microarray', dest='microarray')
	parser.add_argument('-n', '--nanostring', dest='nanostring')
	parser.add_argument('-s', '--sample_conversion', dest='sample_conversion')
	parser.add_argument('-g', '--gene_conversion', dest='gene_conversion')
	parser.add_argument('-t', '--analysis_type', dest='analysis_type',
		choices=['per_gene', 'per_sample'])
	parser.add_argument('-o', '--output_filename', dest='output_filename')
	parsed = parser.parse_args(argv[1:])
	return parsed


def main(argv):
	parsed = parse_args(argv)

	##load expression data
	chip_expr = np.loadtxt(parsed.microarray, dtype=str, delimiter="\t")
	ns_expr = np.loadtxt(parsed.nanostring, dtype=str, delimiter="\t")

	##parse conversion data
	conv_data = np.loadtxt(parsed.sample_conversion, 
		dtype=str, delimiter="\t", usecols=[2,3], skiprows=1)
	chip2ns = {}
	ns2chip = {}
	for i in range(len(conv_data)):
		chip2ns[conv_data[i,0]] = conv_data[i,1]
		ns2chip[conv_data[i,1]] = conv_data[i,0]

	conv_data = np.loadtxt(parsed.gene_conversion, 
		dtype=str, delimiter="\t", usecols=[1,2])
	gene2tc = {}
	for i in range(len(conv_data)):
		tcs = conv_data[i,1].split(",")
		gene2tc[conv_data[i,0]] = [tc for tc in tcs if tc in chip_expr[1:,0]]

	##remove nanostring samples not in microarray data
	valid_ns_sample_indx = [0]
	for j in range(1,ns_expr.shape[1]):
		ns_expr[0,j] = ns_expr[0,j].split(".")[0]
		if ns_expr[0,j] in ns2chip.keys():
			valid_ns_sample_indx.append(j)
	ns_expr = ns_expr[:,valid_ns_sample_indx]

	##remove microarray samples not in nanostring data, and
	##reorder the microarray samples to the corresponding samples in nanostring
	valid_chip_sample_indx = [0]
	for j in range(1,chip_expr.shape[1]):
		chip_expr[0,j] = chip_expr[0,j].split(".")[0]
		if (chip_expr[0,j] in chip2ns.keys()) and (chip2ns[chip_expr[0,j]] in ns_expr[0,1:]):
			valid_chip_sample_indx.append(j)
	chip_expr = chip_expr[:,valid_chip_sample_indx]
	
	ns_sample_indx = [0]
	for i in range(1,ns_expr.shape[1]):
		ns_sample_indx.append(np.where(chip_expr[0,:] == ns2chip[ns_expr[0,i]])[0][0])
	chip_expr = chip_expr[:,ns_sample_indx]

	##take intersection of two gene sets, and average the TCs belonging to the same gene 
	genes = []
	ns_expr_new = []
	chip_expr_new = []
	for i in range(1,ns_expr.shape[0]):
		gene = ns_expr[i,0]
		ns_expr_i = np.array(ns_expr[i,1:], dtype=float)

		if gene in gene2tc.keys():
			tc_indx = []
			for tc in gene2tc[gene]:
				tc_indx.append(np.where(chip_expr[:,0] == tc)[0][0])
			chip_expr_i = np.array(chip_expr[tc_indx,1:], dtype=float)
			chip_expr_i = np.power(2, chip_expr_i) #since microarray expr are in log2
			chip_expr_i = np.mean(chip_expr_i, axis=0)
			# chip_expr_i = np.log2(chip_expr_i)

			genes.append(gene)
			ns_expr_new.append(ns_expr_i)
			chip_expr_new.append(chip_expr_i)
	genes = np.array(genes)
	ns_expr_new = np.array(ns_expr_new)
	chip_expr_new = np.array(chip_expr_new)


	if parsed.analysis_type == 'per_gene':
		##for each gene, compute correlation of all samples
		correlations = []
		for i in range(len(genes)):
			ns_expr_i = ns_expr_new[i,]
			chip_expr_i = chip_expr_new[i,]
			corr = pearsonr(ns_expr_i, chip_expr_i)
			# corr = spearmanr(ns_expr_i, chip_expr_i)
			correlations.append(list(corr))
		correlations = np.array(correlations)

		##sort the correlations
		indx_sorted = np.argsort(correlations[:,0])[::-1]
		correlations = correlations[indx_sorted,]
		genes = genes[indx_sorted]

		out = np.hstack((genes[np.newaxis].T, correlations))
		np.savetxt(parsed.output_filename+".txt", out, fmt="%s", delimiter="\t")


	elif parsed.analysis_type == 'per_sample':
		##for each sample, compute correlation of all genes
		samples = np.array(['_'.join([chip_expr[0,i],ns_expr[0,i]]) 
			for i in range(1,chip_expr.shape[1])])
		
		correlations = []
		for j in range(len(samples)):
			ns_expr_j = ns_expr_new[:,j]
			chip_expr_j = chip_expr_new[:,j]
			corr = pearsonr(ns_expr_j, chip_expr_j)
			# corr = spearmanr(ns_expr_j, chip_expr_j)
			correlations.append(list(corr))

			# print samples[j]
			# if samples[j] == '27562_NS44':
			# 	plt.scatter(ns_expr_j, chip_expr_j)
			# 	plt.show()
			
		correlations = np.array(correlations)

		##sort the correlations
		indx_sorted = np.argsort(correlations[:,0])[::-1]
		correlations = correlations[indx_sorted,]
		samples = samples[indx_sorted]

		out = np.hstack((samples[np.newaxis].T, correlations))
		np.savetxt(parsed.output_filename+".txt", out, fmt="%s", delimiter="\t")

	else:
		sys.exit("Unknown analysis type!")

	##plot histogram
	x = correlations[:,0][np.invert(np.isnan(correlations[:,0]))]
	plt.hist(x, bins=np.arange(-1.,1.,.01), facecolor='blue', alpha=0.75)
	# plt.xlabel('Pearson correlation')
	plt.xlabel('Spearman correlation')
	plt.ylabel('Count')
	plt.show()


if __name__ == "__main__":
    main(sys.argv)
