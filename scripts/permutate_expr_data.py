#!/usr/bin/python
# example:
# python permutate_expr_data.py -d ../data/20160128_project_1930/chipdata_geneset_x_valid_chips.txt -n 1 -o ../data/20160128_project_1930/permutation_experiments

import sys
import os
import argparse
import numpy


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Permutation on gene expression data.')
	parser.add_argument('-d', '--expr_data', dest='expr_data')
	parser.add_argument('-n', '--num_permuted', dest='num_permuted', type=int)
	parser.add_argument('-o', '--output_dir', dest='output_dir')
	parsed = parser.parse_args(argv[1:])
	return parsed


def main(argv):
	parsed = parse_args(argv)
	num_permuted = parsed.num_permuted
	output_dir = parsed.output_dir
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	# read expression data
	expr_data = numpy.loadtxt(parsed.expr_data, dtype=str, delimiter='\t')
	(num_genes, num_samples) = expr_data[1:,1:].shape

	# generate permuted n times
	for i in range(num_permuted):
		# # permute gene expression of each sample
		# permuted_expr_data = numpy.empty((num_genes, num_samples), dtype='|S10')
		# for j in range(num_samples):
		# 	permuted_idx = numpy.random.permutation(num_genes)
		# 	permuted_expr_data[:,j] = expr_data[1:,j+1][permuted_idx]

		# permute gene expression of each sample
		permuted_expr_data = numpy.empty((num_genes, num_samples), dtype='|S10')
		for j in range(num_genes):
			permuted_idx = numpy.random.permutation(num_samples)
			permuted_expr_data[j,:] = expr_data[j+1,1:][permuted_idx]

		# save permuted data
		output_filename = output_dir + '/permuted_expr_data_' + str(i+1) + '.txt'
		output_expr_data = numpy.vstack((expr_data[0],numpy.hstack((expr_data[1:,0][numpy.newaxis].T,permuted_expr_data))))
		numpy.savetxt(output_filename, output_expr_data, fmt='%s', delimiter='\t')


if __name__ == "__main__":
    main(sys.argv)
