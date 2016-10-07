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
	wt_mean = np.mean(expr[:, wt_indices], axis=1)[np.newaxis].T
	return np.log2(expr/wt_mean)


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
	expr = np.power(2, expr)

	## Parse CONTROL sample label
	label = np.empty(len(sample_id), dtype=str)
	for i in range(len(sample_id)):
		if sample_id[i].endswith('.N'):
			label[i] = '0'
		else:
			label[i] = '1'

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
	heatmap(expr_log2_fc[top_indices,:], gene_id[top_indices], sample_id, 'average', 'average', 'euclidean', 'euclidean', 'yellow_black_blue', parsed.dir_output +'top_'+ str(top_rank) +'_'+ parsed.method +'_genes') 
	# heatmap(expr_log2_fc[top_indices,:], limma_genes_for_display, sample_id, 'average', 'average', 'euclidean', 'euclidean', 'yellow_black_blue', parsed.dir_output +'top_'+ str(top_rank) +'_'+ parsed.method +'_genes')


if __name__ == "__main__":
    main(sys.argv)
