#!/usr/bin/python
# example:
# python permutate_de_gene_analysis.py -o ../data/20160128_project_1930/60_limma_genes.txt -d ../data/20160128_project_1930/permutation_experiments/ -n 1 -a ../data/20160128_project_1930/permutation_experiments_results.txt

import sys
import os
import argparse
import numpy
from scipy.stats.mstats import gmean


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Analyze DE genes in permuted expression data.')
	parser.add_argument('-o', '--original_de_genes', dest='original_de_genes')
	parser.add_argument('-d', '--de_genes_dir', dest='de_genes_dir')
	parser.add_argument('-n', '--num_permuted', dest='num_permuted', type=int)
	parser.add_argument('-a', '--output_analysis', dest='output_analysis')
	parsed = parser.parse_args(argv[1:])
	return parsed


def analyze_de_genes(filename):
	data = numpy.loadtxt(filename, delimiter='\t', usecols=range(1,7), skiprows=1)
	(logfc, pval) = (numpy.absolute(data[:,0]), data[:,3])
	return [gmean(logfc), max(logfc), min(logfc), gmean(pval), min(pval), max(pval)]


def main(argv):
	parsed = parse_args(argv)

	# compute geo-mean of absolute logfc, max and min absolute logfc, geo-mean of p-value, max and min p-values  
	results = []
	results.append(analyze_de_genes(parsed.original_de_genes))
	for i in range(parsed.num_permuted):
		results.append(analyze_de_genes(parsed.de_genes_dir +'/60_limma_genes_'+ str(i+1) +'.txt'))

	# write results
	col_names = ['Experiment', 'mean logFC', 'max logFC', 'min logFC', 'mean p-value', 'min p-value', 'max p-value']
	row_names = ['original'] + ['perm_data_'+str(i+1) for i in range(parsed.num_permuted)]
	col_names = numpy.array(col_names, dtype=str)
	row_names = numpy.array(row_names, dtype=str)[numpy.newaxis].T
	results = numpy.array(results, dtype=str)
	results = numpy.vstack((col_names, numpy.hstack((row_names, results))))
	numpy.savetxt(parsed.output_analysis, results, fmt='%s', delimiter='\t')


if __name__ == "__main__":
    main(sys.argv)
