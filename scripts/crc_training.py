#/usr/bin/python

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
	parser = argparse.ArgumentParser(description='Train Random Forest model on EE microarray data.')
	parser.add_argument('-a', '--learning_algorithm', dest='learning_algorithm', default='random_forest', help='options: %s' % learning_algorithms)
	parser.add_argument('-i', '--input_filename', dest='input_filename')
	parser.add_argument('-o', '--output_directory', dest='output_directory')
	parsed = parser.parse_args(argv[1:])
	return parsed

def parse_data(filename, label_col, data_col_start, data_col_end=10000):
	data = np.loadtxt(filename, dtype=str, delimiter='\t')
	gene_id = data[0, data_col_start:data_col_end]
	sample_id = data[1:, 0]
	expr = np.array(data[1:, data_col_start:data_col_end], dtype=np.float32)
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

def main(argv):
	# parse data
	parsed = parse_args(argv)
	if parsed.output_directory != None:
		parsed.output_directory += '/' if (not parsed.output_directory.endswith('/')) else ''
		if (os.path.exists(parsed.output_directory)):
			shutil.rmtree(parsed.output_directory)
		os.makedirs(parsed.output_directory)
	
	[gene_id, sample_id, expr_tr, label_tr] = parse_data(parsed.input_filename, 1, 2)
	
	label_unique= np.unique(label_tr)
	label_count = np.array([len(np.where(label_tr == l)[0]) for l in label_unique])

	print "Training set dimension:", expr_tr.shape[0], "samples x", expr_tr.shape[1], "features"
	print "CRC labels:", label_unique, ", counts:", label_count

	time_start = time.clock()


	##### Random Forest #####
	if parsed.learning_algorithm.lower() ==	'random_forest':
		from sklearn.ensemble import RandomForestClassifier
		
		# train the model
		clf = RandomForestClassifier(n_estimators=1000, n_jobs=1, criterion="entropy", oob_score=True, class_weight=None, verbose=False)
		clf.fit(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + parsed.learning_algorithm.lower() + '_model.pkl')
		
		time_end = time.clock()
		time_elapsed = time_end - time_start

		# sort genes by importance
		num_most_important_gene = 25
		gene_score = clf.feature_importances_
		gene_index = gene_score.argsort()[-num_most_important_gene:][::-1]
		num_most_important_gene = min(num_most_important_gene, len(gene_score))

		# print messages 
		print "Training time elapsed:", time_elapsed, "sec"
		print "Out-of-bag accuracy:", clf.oob_score_


	##### SVM #####
	elif parsed.learning_algorithm.lower() == 'svm':
		from sklearn.svm import SVC

		# train the model
		clf = SVC(C=1.0, kernel='rbf', probability=True, verbose=False)
		clf.fit(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + 'svm_model.pkl')


	##### Neural Network #####
	elif parsed.learning_algorithm.lower() == 'neural_net':
		from sklearn.linear_model import LogisticRegression
		from sklearn.neural_network import BernoulliRBM
		from sklearn.pipeline import Pipeline

		# train the model
		logistic = LogisticRegression(C=10)
		rbm = BernoulliRBM(n_components=256, learning_rate=.001, n_iter=100, verbose=False)
		clf = Pipeline(steps=[('rmb', rbm), ('logistic', logistic)])
		clf.fit(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + parsed.learning_algorithm.lower() + '_model.pkl')


	##### Gradient Boosting #####		
	elif parsed.learning_algorithm.lower() == 'grad_boosting':
		from sklearn.ensemble import GradientBoostingClassifier

		optimal_param = 3
		## cross validation
		# n_folds = 10
		# (expr_tr_cv, label_tr_cv) = generate_cross_validation(expr_tr, label_tr, n_folds=n_folds)
		# param_range = [2,3,4,5,6,7,8,9,10] # choose the param to tune
		# accuracy_lst = []
		# for p in param_range:
		# 	print "Running cross valdiation ... p =", p
		# 	clf = GradientBoostingClassifier(loss='exponential', learning_rate=.0025, n_estimators=1000, max_depth=p, subsample=1.0,verbose=False)
		# 	accuracy_sum = 0
		# 	for i in range(n_folds):
		# 		# internal training and testing
		# 		expr_tr0 = np.vstack(expr_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		label_tr0 = np.hstack(label_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		expr_tr1 = expr_tr_cv[i]
		# 		label_tr1 = label_tr_cv[i]
		# 		clf.fit(expr_tr0, label_tr0)
		# 		label_pred = clf.predict(expr_tr1)
		# 		accuracy_pred = clf.score(expr_tr, label_tr)
		# 		accuracy_sum += accuracy_pred
		# 	accuracy_lst.append(accuracy_sum/float(n_folds))
		# 	print "   Average accuracy:", accuracy_sum/float(n_folds)
		# optimal_param = param_range[np.argmax(accuracy_lst)]
		# print "Optimal param:", optimal_param

		## train the model
		clf = GradientBoostingClassifier(loss='exponential', learning_rate=.0025, n_estimators=1000, max_depth=optimal_param, subsample=1.0,verbose=False)
		clf.fit(expr_tr, label_tr)
		label_pred = clf.predict(expr_tr)
		accuracy_pred = clf.score(expr_tr, label_tr)
		print "Training accuracy:", clf.score(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + parsed.learning_algorithm.lower() + '_model.pkl')

	##### Gaussian Process #####
	elif parsed.learning_algorithm.lower() == 'gauss_process':
		from sklearn.gaussian_process import GaussianProcessClassifier
		from sklearn.gaussian_process.kernels import RBF

		# cross validation
		n_folds = 10
		(expr_tr_cv, label_tr_cv) = generate_cross_validation(expr_tr, label_tr, n_folds=n_folds)
		param_range = [.1,.25,.5,1,2.5,5,10] # choose the param to tune
		accuracy_lst = []
		for p in param_range:
			print "Running cross valdiation ... p =", p
			clf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=p), optimizer="fmin_l_bfgs_b")
			accuracy_sum = 0
			for i in range(n_folds):
				# internal training and testing
				expr_tr0 = np.vstack(expr_tr_cv[np.setdiff1d(range(n_folds),i)])
				label_tr0 = np.hstack(label_tr_cv[np.setdiff1d(range(n_folds),i)])
				expr_tr1 = expr_tr_cv[i]
				label_tr1 = label_tr_cv[i]
				clf.fit(expr_tr0, label_tr0)
				label_pred = clf.predict(expr_tr1)
				accuracy_pred = len([label_pred[i] for i in range(len(label_pred)) if (label_pred[i] == label_tr1[i])]) / float(len(label_pred))
				accuracy_sum += accuracy_pred
			accuracy_lst.append(accuracy_sum/float(n_folds))
			print "   Average accuracy:", accuracy_sum/float(n_folds)
		optimal_param = param_range[np.argmax(accuracy_lst)]
		print "Optimal param:", optimal_param

		# train the model
		clf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=1.0), optimizer="fmin_l_bfgs_b")
		clf.fit(expr_tr, label_tr)
		label_pred = clf.predict(expr_tr)
		print "Training accuracy:", clf.score(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + parsed.learning_algorithm.lower() + '_model.pkl')

	else:
		sys.exit('Improper learning algorithm option given.')

	# print timing messages
	time_end = time.clock()
	time_elapsed = time_end - time_start
	print "Training time elapsed:", time_elapsed, "sec"


if __name__ == "__main__":
    main(sys.argv)

