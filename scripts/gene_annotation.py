#!/usr/bin/python
# example:
#

import sys
import argparse
import numpy


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Look up gene annotations.')
	parser.add_argument('-t', '--top_genes', dest='top_genes')
	parser.add_argument('-d', '--expr_data', dest='expr_data')
	parser.add_argument('-n', '--num_samples', dest='num_samples', type=int)
	parser.add_argument('-p', '--pval', dest='pval', type=float)
	parser.add_argument('-o', '--output_annotation', dest='output_annotation')
	parsed = parser.parse_args(argv[1:])
	return parsed


def main(argv):
	parsed = parse_args(argv)

	# load data
	num_samples = parsed.num_samples
	top_genes = numpy.loadtxt(parsed.top_genes, dtype=str, skiprows=1)
	expr_data = numpy.loadtxt(parsed.expr_data, dtype=str, delimiter='\t')

	# annotations: 'Probe Set ID' 'Gene Accession' 'Gene Symbol' 'Gene Description' 'mRNA Accession' 'mRna - Description'
	col_idx = [0] + range(num_samples+1, num_samples+4) + [num_samples+5] + [num_samples+7]
	annot = expr_data[0][col_idx]
	
	if parsed.pval == None:
		print "p-value is None."
		for i in range(len(top_genes)):
			row_idx = numpy.where(expr_data[:,0] == top_genes[i,0])[0][0]
			annot = numpy.vstack((annot, expr_data[row_idx][col_idx]))
	else:
		print "p-value is not None."
		pval = parsed.pval	
		for i in range(len(top_genes)):
			if float(top_genes[i,4]) < pval:
				row_idx = numpy.where(expr_data[:,0] == top_genes[i,0])[0][0]
				annot = numpy.vstack((annot, expr_data[row_idx][col_idx]))
			else:
				break

	# save data
	numpy.savetxt(parsed.output_annotation, annot, fmt='%s', delimiter='\t')


if __name__ == "__main__":
    main(sys.argv)
