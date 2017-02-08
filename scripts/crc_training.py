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
	parser.add_argument('-a', '--learning_algorithm', dest='learning_algorithm', default='random_forest', help='options: %s' % learning_algorithms)
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


	##### Random Forest #####
	if parsed.learning_algorithm.lower() ==	'random_forest':
		from sklearn.ensemble import RandomForestClassifier
		
		# train the model
		clf = RandomForestClassifier(n_estimators=1000, n_jobs=5, criterion="entropy", oob_score=True, class_weight="balanced", verbose=False)
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

		# ## cross validation
		# n_folds = 10
		# (expr_tr_cv, label_tr_cv) = generate_cross_validation(expr_tr, label_tr, n_folds=n_folds)
		# param_range = range(1,8) ##max depth
		# # param_range = [.0005, .001, .0015, .002, .0025, .005, .0075, .01] ##learning rate
		# # param_range = np.arange(.1,1,.1) ##stochastic process
		# accuracy_lst = []
		# for p in param_range:
		# 	print "Running cross valdiation ... p =", p
		# 	clf = GradientBoostingClassifier(learning_rate=.0025, n_estimators=1000, max_depth=p, subsample=1.0,verbose=False)
		# 	accuracy_sum = 0
		# 	for i in range(n_folds):
		# 		# internal training and testing
		# 		expr_tr0 = np.vstack(expr_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		label_tr0 = np.hstack(label_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		expr_tr1 = expr_tr_cv[i]
		# 		label_tr1 = label_tr_cv[i]
		# 		clf.fit(expr_tr0, label_tr0)
		# 		label_pred = clf.predict(expr_tr1)
		# 		# accuracy_pred = clf.score(expr_tr, label_tr)
		# 		[foo, foo, accuracy_pred] = calculate_confusion_matrix(label_tr1, label_pred)
		# 		accuracy_sum += accuracy_pred
		# 	accuracy_lst.append(accuracy_sum/float(n_folds))
		# 	print "   Average accuracy:", accuracy_sum/float(n_folds)
		# optimal_param = param_range[np.argmax(accuracy_lst)]
		# print "Optimal param:", optimal_param

		## train the model
		# clf = GradientBoostingClassifier(learning_rate=.0025, n_estimators=1000, max_depth=optimal_param, subsample=1.0,verbose=False)
		clf = GradientBoostingClassifier(learning_rate=.01, n_estimators=1000, max_depth=optimal_param, subsample=0.8, verbose=False)

		clf.fit(expr_tr, label_tr)
		label_pred = clf.predict(expr_tr)
		accuracy_pred = clf.score(expr_tr, label_tr)
		print "Training accuracy:", clf.score(expr_tr, label_tr)
		if parsed.output_directory != None:
			joblib.dump(clf, parsed.output_directory + parsed.learning_algorithm.lower() + '_model.pkl')

		# ## calculate score for ML prediction
		# score_ml = .8*clf.predict_proba(expr_tr)[:,0]
		# ## add weight to outlier predictors
		# score_predictors = np.zeros(len(sample_id))
		# (tc_predictors, tc_predictors_expr) = get_predictor_expr(parsed.outlier_predictors, expr_tr_full, gene_id_full)
		# print "Predictors added:", tc_predictors
		# weight_predictors = [.2/len(tc_predictors)]*len(tc_predictors)
		# (tc_predictors_normal_expr_median, tc_predictors_normal_expr_sd, tc_predictors_first_pct, tc_predictors_third_pct) = parse_predictor_stats(tc_predictors_expr)
		# # thlds = tc_predictors_normal_expr_median + 2*tc_predictors_normal_expr_sd
		# thlds = 1*(tc_predictors_third_pct - tc_predictors_first_pct) + tc_predictors_third_pct
		# for j in range(len(tc_predictors)):
		# 	indx = np.where(tc_predictors_expr[:,j] > thlds[j])[0]
		# 	for i in indx:
		# 		score_predictors[i] += weight_predictors[j]
		# 		# score_predictors[i] = .2
		# ## save normal statistis 
		# out_normal_stats = [["TC_id", "median", "sd", "1st_pct", "3rd_pct"]]
		# for j in range(len(tc_predictors)):
		# 	out_normal_stats.append([tc_predictors[j], tc_predictors_normal_expr_median[j], tc_predictors_normal_expr_sd[j], tc_predictors_first_pct[j], tc_predictors_third_pct[j]])
		# np.savetxt(parsed.normal_stats, np.array(out_normal_stats, dtype=str), fmt="%s", delimiter="\t")
		# ## final prediction
		# print("\t".join(["sample_id", "true_label", "predicted_label", "score_ML", "score_predictors", "final_score", "score_change?"]))
		# score_final = []
		# for i in range(len(sample_id)):
		# 	score = (score_ml[i]+score_predictors[i])/(1+score_predictors[i])
		# 	score_final.append(score)
		# 	predicted_label = "C" if score > .5 else "N"
		# 	print("\t".join([sample_id[i], label_tr[i], predicted_label, str(score_ml[i]), str(score_predictors[i]), str(score), str(score != score_ml[i])] ))
		


	##### AdaBoost #####
	elif parsed.learning_algorithm.lower() == "adaboost":
		from sklearn.ensemble import AdaBoostClassifier
		from sklearn.tree import DecisionTreeClassifier
		clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), n_estimators=1000, learning_rate=.0025)
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

		optimal_param = 1.0
		## cross validation
		# n_folds = 10
		# (expr_tr_cv, label_tr_cv) = generate_cross_validation(expr_tr, label_tr, n_folds=n_folds)
		# param_range = [.1,.25,.5,1,2.5,5,10] # choose the param to tune
		# accuracy_lst = []
		# for p in param_range:
		# 	print "Running cross valdiation ... p =", p
		# 	clf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=p), optimizer="fmin_l_bfgs_b")
		# 	accuracy_sum = 0
		# 	for i in range(n_folds):
		# 		# internal training and testing
		# 		expr_tr0 = np.vstack(expr_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		label_tr0 = np.hstack(label_tr_cv[np.setdiff1d(range(n_folds),i)])
		# 		expr_tr1 = expr_tr_cv[i]
		# 		label_tr1 = label_tr_cv[i]
		# 		clf.fit(expr_tr0, label_tr0)
		# 		label_pred = clf.predict(expr_tr1)
		# 		accuracy_pred = len([label_pred[i] for i in range(len(label_pred)) if (label_pred[i] == label_tr1[i])]) / float(len(label_pred))
		# 		accuracy_sum += accuracy_pred
		# 	accuracy_lst.append(accuracy_sum/float(n_folds))
		# 	print "   Average accuracy:", accuracy_sum/float(n_folds)
		# optimal_param = param_range[np.argmax(accuracy_lst)]
		# print "Optimal param:", optimal_param

		# train the model
		clf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=optimal_param), optimizer="fmin_l_bfgs_b")
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

