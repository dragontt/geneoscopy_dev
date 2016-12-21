#!/usr/bin/python
import sys
import argparse
import numpy as np
import random

def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-i', '--input_expr', dest='input_expr')
	parser.add_argument('-v', '--valid_chips', dest='valid_chips')
	parser.add_argument('-f', '--n_folds', dest='n_folds', type=int, default=10)
	parser.add_argument('-g', '--group', dest='group')
	parser.add_argument('-o', '--output_directory', dest='output_directory')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	valid = np.loadtxt(parsed.valid_chips, dtype=str, delimiter="\t")
	expr = np.loadtxt(parsed.input_expr, dtype=str, delimiter="\t")
	labels = ['header']
	for i in range(1,expr.shape[1]):
		labels.append(expr[0,i].split(".")[1])
	labels = np.array(labels)

	if parsed.group == 'N_vs_C':
		indx_N = np.where(labels == "N")[0]
		indx_C = np.where(labels == "C")[0]
		random.shuffle(indx_N)
		random.shuffle(indx_C)
		indx_N_rand_chunks = np.array_split(indx_N, parsed.n_folds)
		indx_C_rand_chunks = np.array_split(indx_C, parsed.n_folds)
		for i in range(parsed.n_folds):
			indx_combined = [0]
			for j in range(parsed.n_folds):
				if j != i:
					indx_combined += list(indx_N_rand_chunks[j]) 
					indx_combined += list(indx_C_rand_chunks[j])
			rand_expr = expr[:,indx_combined]
			del indx_combined[0]
			indx_combined = [k-1 for k in indx_combined]
			rand_valid = valid[indx_combined]
			np.savetxt(parsed.output_directory+'/chipdata_random_'+str(i+1)+'.txt',
				rand_expr, fmt='%s', delimiter='\t')
			np.savetxt(parsed.output_directory+'/valid_chips_random_'+str(i+1)+'.txt',
				rand_valid, fmt='%s', delimiter='\t')

	elif parsed.group == 'N_vs_P_vs_C':
		indx_N = np.where(labels == "N")[0]
		indx_P = np.where(labels == "P")[0]
		indx_C = np.where(labels == "C")[0]
		random.shuffle(indx_N)
		random.shuffle(indx_P)
		random.shuffle(indx_C)
		indx_N_rand_chunks = np.array_split(indx_N, parsed.n_folds)
		indx_P_rand_chunks = np.array_split(indx_P, parsed.n_folds)
		indx_C_rand_chunks = np.array_split(indx_C, parsed.n_folds)
		for i in range(parsed.n_folds):
			indx_combined = [0]
			for j in range(parsed.n_folds):
				if j != i:
					indx_combined += list(indx_N_rand_chunks[j]) 
					indx_combined += list(indx_P_rand_chunks[j])
					indx_combined += list(indx_C_rand_chunks[j])
			rand_expr = expr[:,indx_combined]
			del indx_combined[0]
			indx_combined = [k-1 for k in indx_combined]
			rand_valid = valid[indx_combined]
			np.savetxt(parsed.output_directory+'/chipdata_random_'+str(i+1)+'.txt',
				rand_expr, fmt='%s', delimiter='\t')
			np.savetxt(parsed.output_directory+'/valid_chips_random_'+str(i+1)+'.txt',
				rand_valid, fmt='%s', delimiter='\t')

	else:
		sys.exit("Group unavalilable!")


if __name__ == "__main__":
    main(sys.argv)