#!/usr/bin/python
import numpy as np
from scipy.stats import mannwhitneyu

data = np.loadtxt('/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_lit.txt', dtype=str, delimiter='\t')
samples = data[0,6:]
groups = data[1,6:]
genes = data[2:,1]
expr = np.array(data[2:,6:], dtype=float)

normal_indx = np.where(groups == "Normal")[0]
polyp_indx = np.where(groups == "Polyp")[0]
cancer_indx = np.where(groups == "Cancer")[0]

for contrast_group_indx in [np.append(polyp_indx,cancer_indx)]:
# for contrast_group_indx in [polyp_indx, cancer_indx, np.append(polyp_indx,cancer_indx)]:
	utest_results = []
	for i in range(len(genes)):
		control = expr[i,normal_indx]
		contrast = expr[i,contrast_group_indx]
		utest_stats = mannwhitneyu(control, contrast)
		fc = np.mean(contrast)/np.mean(control)
		utest_results.append([np.log2(fc), utest_stats[1]])
	utest_results = np.array(utest_results)

	indx_sorted = np.argsort(utest_results[:,1])
	# indx_sorted = np.argsort(np.abs(utest_results[:,0]))[::-1]
	out = np.hstack(( genes[indx_sorted][np.newaxis].T, utest_results[indx_sorted,:]))
	np.savetxt('/Users/KANG/geneoscopy_dev/data/20170113_nanostring_project_18731/POP_48_samples_011817_PosNormData_lit_DE_analysis.txt', out, fmt="%s", delimiter="\t")

	# print "gene", "log2FC", "p-val"

