#!/usr/bin/python
import sys
import os
import argparse
import shutil
import numpy as np
import time
from sklearn.externals import joblib
import random

learning_algorithms = ['random_forest', 'svm', 'neural_net', 'grad_boosting']

def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-i', '--input_expr', dest='input_expr')
	parser.add_argument('-f', '--input_expr_full', dest='input_expr_full')
	parser.add_argument('-p', '--outlier_predictors', dest='outlier_predictors')
	parser.add_argument('-s', '--normal_stats', dest='normal_stats')
	parser.add_argument('-o', '--output_directory', dest='output_directory')
	parsed = parser.parse_args(argv[1:])
	return parsed


def parse_data(filename, label_col, data_col_start):
	data = np.loadtxt(filename, dtype=str, delimiter='\t')
	gene_id = data[0, data_col_start:]
	sample_id = data[1:, 0]
	expr = np.array(data[1:, data_col_start:], dtype=np.float32)
	label = data[1:, label_col]
	return [gene_id, sample_id, expr, label]


def generate_cross_validation(expr_tr, label_tr, n_folds=10):
	# make sure at least one Normal, at least one CRC in each fold
	indx_N = np.where(label_tr == "N")[0]
	indx_C = np.setdiff1d(range(len(label_tr)), indx_N)
	random.shuffle(indx_N)
	random.shuffle(indx_C)
	indx_N_rand_chunks = np.array_split(indx_N, n_folds)
	indx_C_rand_chunks = np.array_split(indx_C, n_folds)
	expr_tr_rand_chuks = []
	label_tr_rand_chuks = []
	for i in range(len(indx_N_rand_chunks)):
		indx_combined = list(indx_N_rand_chunks[i]) + list(indx_C_rand_chunks[i])
		expr_tr_rand_chuks.append(expr_tr[indx_combined,:])
		label_tr_rand_chuks.append(label_tr[indx_combined])
	return (np.array(expr_tr_rand_chuks), np.array(label_tr_rand_chuks))


def get_predictor_expr(filename, expr, gene_id):
	f = open(filename, "r")
	lines = f.readlines()
	tc_predictors_indx = []
	tc_predictors = []
	for line in lines:
		tmp_arr = line.strip().split("\t")
		if len(tmp_arr) > 1:
			for tmp_tc in tmp_arr[1].split(","):
				if tmp_tc in gene_id:
					tmp_indx = np.where(gene_id == tmp_tc)[0][0]
					tc_predictors_indx.append(tmp_indx)
					tc_predictors.append(tmp_tc)
	return (tc_predictors, expr[:, tc_predictors_indx])


def parse_predictor_stats(expr):
	return (np.median(expr, axis=0), np.median(expr, axis=0), np.percentile(expr, 25, axis=0),np.percentile(expr, 75, axis=0))


def boostrap_label_group(sample_id, expr, labels):
	## boostrap label groups that have smaller number of samples to match the largest
	unique_labels = np.unique(labels)
	label_indx_dict = {}
	largest_sample_num = 0
	for ul in unique_labels:
		tmp_arr = np.where(labels == ul)[0]
		label_indx_dict[ul] = tmp_arr
		largest_sample_num = len(tmp_arr) if len(tmp_arr)>largest_sample_num else largest_sample_num
	for ul in unique_labels:
		tmp_arr = label_indx_dict[ul]
		if len(tmp_arr) != largest_sample_num:
			new_arr = np.random.choice(tmp_arr, largest_sample_num-len(tmp_arr), replace=True)
			sample_id = np.append(sample_id, sample_id[new_arr])
			expr = np.vstack((expr, expr[new_arr,]))
			labels = np.append(labels, labels[new_arr])
	return (sample_id, expr, labels)


def subsample_label_group(sample_id, expr, labels):
	## subsample label groups that have larger number of samples to match the smallest
	unique_labels = np.unique(labels)
	label_indx_dict = {}
	smallest_sample_num = 10000000000
	for ul in unique_labels:
		tmp_arr = np.where(labels == ul)[0]
		label_indx_dict[ul] = tmp_arr
		smallest_sample_num = len(tmp_arr) if len(tmp_arr)<smallest_sample_num else smallest_sample_num
	indx_to_remove = np.array([])
	for ul in unique_labels:
		tmp_arr = label_indx_dict[ul]
		if len(tmp_arr) != smallest_sample_num:
			new_arr = np.random.choice(tmp_arr, len(tmp_arr)-smallest_sample_num, replace=False)
			indx_to_remove = np.append(indx_to_remove, new_arr)
	sample_id = np.delete(sample_id, indx_to_remove)
	expr = np.delete(expr, indx_to_remove, axis=0)
	labels = np.delete(labels, indx_to_remove)
	return (sample_id, expr, labels)


def calculate_confusion_matrix(label_te, label_pred):
	pred_pos = np.append(np.where(label_pred == "C")[0], np.where(label_pred == "P")[0])
	pred_neg = np.where(label_pred == "N")[0]
	te_pos = np.append(np.where(label_te == "C")[0], np.where(label_te == "P")[0])
	te_neg = np.where(label_te == "N")[0]
	tps = len(np.intersect1d(pred_pos, te_pos)) 
	fps = len(np.intersect1d(pred_pos, te_neg)) 
	fns = len(np.intersect1d(pred_neg, te_pos))
	tns = len(np.intersect1d(pred_neg, te_neg))
	print tps, fps, fns, tns
	sens = tps/float(len(te_pos))
	spec = tns/float(len(te_neg))
	accu = (tps+tns)/float(len(label_te))
	return [sens, spec, accu]


def main(argv):
	# parse data
	parsed = parse_args(argv)
	if parsed.output_directory != None:
		parsed.output_directory += '/' if (not parsed.output_directory.endswith('/')) else ''
		if (os.path.exists(parsed.output_directory)):
			shutil.rmtree(parsed.output_directory)
		os.makedirs(parsed.output_directory)
	
	[gene_id, sample_id, expr_tr, label_tr] = parse_data(parsed.input_expr, 1, 2)
	[gene_id_full, foo, expr_tr_full, foo] = parse_data(parsed.input_expr_full, 1, 2)

	## boostrap the label groups with smaller samples
	# [foo, expr_tr_full, foo] = boostrap_label_group(sample_id, expr_tr_full, label_tr)
	# [sample_id, expr_tr, label_tr] = boostrap_label_group(sample_id, expr_tr, label_tr)
	## subsample the label groups with larger samples
	# [foo, expr_tr_full, foo] = subsample_label_group(sample_id, expr_tr_full, label_tr)
	# [sample_id, expr_tr, label_tr] = subsample_label_group(sample_id, expr_tr, label_tr)
	
	label_unique= np.unique(label_tr)
	label_count = np.array([len(np.where(label_tr == l)[0]) for l in label_unique])

	print "Training set dimension:", expr_tr.shape[0], "samples x", expr_tr.shape[1], "features"
	print "CRC labels:", label_unique, ", counts:", label_count

	time_start = time.clock()

	##### ENSEMBLE MODEL #####
	from sklearn.ensemble import VotingClassifier
	from sklearn.ensemble import RandomForestClassifier
	from sklearn.svm import SVC
	from sklearn.ensemble import GradientBoostingClassifier
	from sklearn.ensemble import AdaBoostClassifier
	from sklearn.tree import DecisionTreeClassifier
	from sklearn.gaussian_process import GaussianProcessClassifier
	from sklearn.gaussian_process.kernels import RBF
	from sklearn.model_selection import cross_val_score

	clfrf = RandomForestClassifier(n_estimators=1000)
	clfsvm = SVC(C=1.6, 
						kernel='rbf',
						probability=True)
	clfgb = GradientBoostingClassifier(learning_rate=.0025, 
						max_depth=3,
						subsample=.8,
						n_estimators=1000)
	clfab = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3))
	clfgp = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=1.0), 
						optimizer="fmin_l_bfgs_b")

	## define model
	eclf = VotingClassifier(estimators=[('rf',clfrf), ('svm',clfsvm), ('gb',clfgb), ('ab',clfab), ('gp',clfgp)], voting="soft", weights=[3,5,1,1,3])

	## cross validation
	# for clf, label in zip([clfrf, clfsvm, clfgb, clfab, clfgp, eclf], 
	# 	['Random forest', 'SVM', 'Grad boost', 'Adabost', 'Gauss proc', 'Ensemble']):
	# 	scores = cross_val_score(clf, expr_tr, label_tr, cv=10, scoring='accuracy')
	# 	print("Accuracy: %0.5f (+/- %0.5f) [%s]" % (scores.mean(), scores.std(), label))
	
	## fit model
	eclf = eclf.fit(expr_tr, label_tr)

	## save the model
	if parsed.output_directory != None:
		joblib.dump(eclf, parsed.output_directory + '/ensemble_model.pkl')


	## print timing messages
	time_end = time.clock()
	time_elapsed = time_end - time_start
	print "Training time elapsed:", time_elapsed, "sec"


if __name__ == "__main__":
    main(sys.argv)

