#! /usr/bin/python
import sys
import argparse
from hierarchical_clustering import heatmap
import matplotlib.lines as plin
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import numpy as np
from scipy.stats.mstats import gmean


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Compute and plot hierachical clustering.')
	parser.add_argument('-m', '--method', dest='method', default='limma')
	parser.add_argument('-i', '--input_expr', dest='input_expr')
	parser.add_argument('-l', '--limma_genes', dest='limma_genes')
	parser.add_argument('-p', '--limma_pval', dest='limma_pval', type=float)
	parser.add_argument('-r', '--top_rank', dest='top_rank', type=int, default=50)
	parser.add_argument('-o', '--dir_output', dest='dir_output')
	parsed = parser.parse_args(argv[1:])
	return parsed


def compute_log2_fc(expr, label):
	wt_indices = np.where(label=='0')[0]
	# wt_mean = np.mean(expr[:, wt_indices], axis=1)[np.newaxis].T
	wt_median = np.median(expr[:, wt_indices], axis=1)[np.newaxis].T
	return np.log2(expr/wt_median)


def compute_fc(expr, label):
	wt_indices = np.where(label=='0')[0]
	# wt_mean = np.mean(expr[:, wt_indices], axis=1)[np.newaxis].T
	wt_median = np.median(expr[:, wt_indices], axis=1)[np.newaxis].T
	return expr/wt_median


def rank_de_genes(expr, label):
	cancer_indices = np.where(label=='1')[0]
	wt_indices = np.where(label=='0')[0]
	cancer_mean = np.mean(expr[:,cancer_indices], axis=1)
	wt_mean = np.mean(expr[:,wt_indices], axis=1)
	return np.log2(cancer_mean/wt_mean)


def fld_de_genes(expr, label):
	cancer_indices = np.where(label=='1')[0]
	wt_indices = np.where(label=='0')[0]
	cancer_mean = np.mean(expr[:,cancer_indices], axis=1)
	wt_mean = np.mean(expr[:,wt_indices], axis=1)
	cancer_sd = np.std(expr[:,cancer_indices], axis=1)
	wt_sd = np.std(expr[:,wt_indices], axis=1)
	fld = np.power(cancer_mean-wt_mean,2)/(np.power(cancer_sd,2)+np.power(wt_sd,2))
	return (fld, ~np.isnan(fld))


def main(argv):
	parsed = parse_args(argv)
	parsed.dir_output += '/'

	## Parse expr data
	raw_data = np.loadtxt(parsed.input_expr, dtype=str, delimiter='\t')
	gene_id = raw_data[1:, 0]
	sample_id = raw_data[0, 1:]

	## Transform from log2 to regular expr values
	expr = np.array(raw_data[1:, 1:], dtype=float)
	# expr = np.array(raw_data[1:, 2:], dtype=float)
	expr = np.power(2, expr)

	## Parse CONTROL sample label
	label = np.empty(len(sample_id), dtype=str)
	label_N = []
	label_P = []
	label_C = []
	for i in range(len(sample_id)):
		if sample_id[i].endswith('.N'):
			label[i] = '0'
			label_N.append(i)
		elif sample_id[i].endswith('.P'):
			label[i] = '1'
			label_P.append(i)
		elif sample_id[i].endswith('.C'):
			label[i] = '1'
			label_C.append(i)

	## Differnetial expression analyssi
	if parsed.method.lower() == 'limma': # Limma de analysis  
		limma_genes = np.loadtxt(parsed.limma_genes, dtype=str, skiprows=1, usecols=[0])
		limma_genes_for_display = np.loadtxt(parsed.limma_genes, dtype=str, skiprows=1, usecols=[4])
		if parsed.limma_pval == None:
			top_rank = len(limma_genes)
		else:
			pvals = np.loadtxt(parsed.limma_genes, skiprows=1, usecols=[4])
			top_rank = np.where(pvals < parsed.limma_pval)[0][-1] + 1
		top_indices = np.zeros(top_rank, dtype=int)
		for i in range(top_rank):
			top_indices[i] = np.where(gene_id == limma_genes[i])[0][0]
		expr_log2_fc = compute_log2_fc(expr, label)
		# expr_fc = compute_fc(expr, label)

	elif parsed.method.lower() == 'fld': # Fisher's linear descriminant
		top_rank = parsed.top_rank
		# rank gene based on fld  
		(de_genes, valid_boolean) = fld_de_genes(expr, label)
		# clean up the nan FLD scores
		de_genes = de_genes[valid_boolean]
		gene_id = gene_id[valid_boolean]
		expr = expr[valid_boolean]
		top_indices = np.argsort(de_genes)[::-1][:top_rank]
		np.savetxt(parsed.dir_output +'top_'+ str(top_rank) +'_'+ parsed.method +'_genes.txt', np.transpose(np.vstack((gene_id[top_indices], de_genes[top_indices]))), fmt='%s', delimiter='\t')
		expr_log2_fc = compute_log2_fc(expr, label)

	elif parsed.method.lower() == 'straight_de':
		top_rank = parsed.top_rank
		# rank gene based on absolute log2 fc 
		de_genes = rank_de_genes(expr, label)
		top_indices = np.argsort(np.absolute(de_genes))[::-1][:top_rank]
		np.savetxt(parsed.dir_output +'top_'+ str(top_rank) +'_'+ parsed.method +'_genes.txt', np.transpose(np.vstack((gene_id[top_indices], de_genes[top_indices]))), fmt='%s', delimiter='\t')
		expr_log2_fc = compute_log2_fc(expr, label)

	else: 
		sys.exit('No DE analysis method given.')

	## Analyze individual DE fold change 
	val_thld = 1.2 #.8
	expr_to_plot = expr_log2_fc[top_indices,:]
	# expr_to_plot = expr_fc[top_indices,:]
	expr_to_plot[expr_to_plot > val_thld] = val_thld
	expr_to_plot[expr_to_plot < -val_thld] = -val_thld

	label_C_subindx = [10, 24, 71, 87, 8, 13, 9, 0, 1, 73, 80, 34, 38, 70, 78, 20, 88, 28, 63, 18, 79, 43, 4, 26, 52, 68, 31, 67, 86, 81, 33, 84, 64, 37, 41, 83, 74, 93, 23, 75, 2, 92, 82, 14, 51, 57, 7, 17, 5, 30, 77, 15, 21, 47, 48, 72, 22, 89, 66, 69, 39, 90, 45, 32, 61, 19, 65, 49, 54, 60, 27, 40, 42, 29, 85, 91, 46, 44, 25, 58, 35, 36, 53, 56, 62, 59, 6, 76, 50, 55, 16, 3, 11, 12]
	label_P_subindx = [14, 68, 77, 34, 78, 60, 71, 69, 65, 4, 52, 10, 1, 67, 75, 0, 53, 66, 79, 26, 73, 12, 28, 32, 46, 9, 15, 19, 23, 25, 11, 63, 17, 45, 44, 5, 51, 18, 56, 6, 58, 70, 41, 72, 59, 62, 64, 36, 48, 61, 24, 76, 50, 74, 27, 30, 39, 31, 55, 13, 3, 7, 81, 16, 33, 40, 57, 21, 49, 54, 22, 47, 80, 2, 38, 42, 29, 35, 37, 43, 8, 20]	
	label_N_subindx = [37, 64, 21, 40, 6, 31, 19, 23, 27, 17, 74, 69, 26, 36, 29, 84, 46, 56, 48, 61, 49, 68, 13, 60, 0, 14, 20, 7, 18, 53, 54, 58, 33, 85, 65, 52, 88, 86, 76, 75, 81, 83, 11, 63, 80, 10, 66, 72, 28, 67, 35, 44, 4, 38, 73, 32, 59, 78, 57, 79, 39, 41, 30, 45, 47, 22, 24, 12, 15, 71, 55, 70, 82, 9, 51, 87, 43, 77, 5, 50, 3, 16, 25, 34, 42, 62, 8, 1, 2]
	indx_samples = np.append(np.append(
					np.array(label_C, dtype=int)[label_C_subindx],
					np.array(label_P, dtype=int)[label_P_subindx]),
					np.array(label_N, dtype=int)[label_N_subindx])

	# indx_samples = label_C + label_P + label_N
	# indx_samples = label_C

	expr_to_plot = expr_to_plot[:, indx_samples]
	sample_id = sample_id[indx_samples]

	# heatmap(expr_to_plot, gene_id[top_indices], sample_id, 'average', 'average', 'euclidean', 'euclidean', 'yellow_black_blue', parsed.dir_output +'top_'+ str(top_rank) +'_'+ parsed.method +'_genes') 
	heatmap(expr_to_plot, gene_id[top_indices], sample_id, 'average', None, 'euclidean', None, 'yellow_black_blue', parsed.dir_output +'new.top_'+ str(top_rank) +'_'+ parsed.method +'_genes') 


if __name__ == "__main__":
    main(sys.argv)
