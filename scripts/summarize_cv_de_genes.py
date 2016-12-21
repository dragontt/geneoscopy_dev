#!/usr/bin/python
import sys
import argparse
import numpy as np
# import scipy.stats as ss
import operator

def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-i', '--cv_directory', dest='cv_directory')
	parser.add_argument('-f', '--n_folds', dest='n_folds', type=int, default=10)
	parser.add_argument('-t', '--n_top_genes', dest='n_top_genes', type=int)
	parser.add_argument('-o', '--output_file', dest='output_file')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	rankings = {}
	for i in range(parsed.n_folds):
		cv_de = np.loadtxt(parsed.cv_directory+'/top_de_genes_'+str(i+1)+'.txt', dtype=str, skiprows=1, usecols=[0])
		if i == 0:
			for j in range(len(cv_de)):
				rankings[cv_de[j]] = j
		else:
			for j in range(len(cv_de)):
				rankings[cv_de[j]] += j
	
	out = np.array(sorted(rankings.items(), key=operator.itemgetter(1)))
	header = np.array(["gene", "sum_ranking"])[np.newaxis]
	out = np.vstack((header, out[:parsed.n_top_genes]))
	np.savetxt(parsed.output_file, out, fmt="%s", delimiter="\t")

if __name__ == "__main__":
    main(sys.argv)