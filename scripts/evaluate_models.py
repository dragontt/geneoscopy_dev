#!/usr/bin/python
import numpy as np
import json
import os
import time
import sys

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.svm import NuSVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB


def split_expr(expr, label_split):
	label = expr[1:,1]
	expr = np.array(expr[1:,2:], dtype=np.float32)
	indx_tr = []
	for k in label_split['training'].keys():
		indx_tr += label_split['training'][k]
	indx_te = []
	for k in label_split['testing'].keys():
		indx_te += label_split['testing'][k]
	return [expr[indx_tr,:], label[indx_tr], expr[indx_te,:], label[indx_te]]


def write_valid_chips(label_indx_dict, samples, labels, filename):
	out = []
	for i in label_indx_dict['N']:
		out.append([samples[i]+'.'+labels[i], '0', '0'])
	for i in label_indx_dict['P']:
		out.append([samples[i]+'.'+labels[i], '1', '0'])
	for i in label_indx_dict['C']:
		out.append([samples[i]+'.'+labels[i], '2', '0'])
	np.savetxt(filename, np.array(out), fmt='%s', delimiter='\t')


def load_de_genes(filename):
	return np.loadtxt(filename, dtype=str, skiprows=1, usecols=[0])


def filter_features(expr, features):
	indx = np.array([0,1])
	for f in features:
		indx = np.append(indx, np.where(expr[0,:]==f)[0][0])
	return expr[:,indx]


def calculate_confusion_matrix(label_te, label_pred):
	pred_pos = np.append(np.where(label_pred == "C")[0], np.where(label_pred == "P")[0])
	pred_neg = np.where(label_pred == "N")[0]
	te_pos = np.append(np.where(label_te == "C")[0], np.where(label_te == "P")[0])
	te_neg = np.where(label_te == "N")[0]
	tps = len(np.intersect1d(pred_pos, te_pos)) 
	fps = len(np.intersect1d(pred_pos, te_neg)) 
	fns = len(np.intersect1d(pred_neg, te_pos))
	tns = len(np.intersect1d(pred_neg, te_neg))
	# print tps, fps, fns, tns
	sens = tps/float(len(te_pos))
	spec = tns/float(len(te_neg))
	accu = (tps+tns)/float(len(label_te))
	return [sens, spec, accu]


def calculate_poly_detection(label_te, label_pred):
	pred_pos = np.append(np.where(label_pred == "P")[0], np.where(label_pred == "C")[0])
	pred_neg = np.where(label_pred == "N")[0]
	te_pos = np.where(label_te == "P")[0]
	te_neg = np.where(label_te == "N")[0]
	tps = len(np.intersect1d(pred_pos, te_pos)) 
	fps = len(np.intersect1d(pred_pos, te_neg)) 
	fns = len(np.intersect1d(pred_neg, te_pos))
	tns = len(np.intersect1d(pred_neg, te_neg))
	sens = tps/float(len(te_pos))
	spec = tns/float(len(te_neg))
	return [sens, spec]


def calculate_cancer_detection(label_te, label_pred):
	pred_pos = np.append(np.where(label_pred == "P")[0], np.where(label_pred == "C")[0])
	pred_neg = np.where(label_pred == "N")[0]
	te_pos = np.where(label_te == "C")[0]
	te_neg = np.where(label_te == "N")[0]
	tps = len(np.intersect1d(pred_pos, te_pos)) 
	fps = len(np.intersect1d(pred_pos, te_neg)) 
	fns = len(np.intersect1d(pred_neg, te_pos))
	tns = len(np.intersect1d(pred_neg, te_neg))
	sens = tps/float(len(te_pos))
	spec = tns/float(len(te_neg))
	return [sens, spec]


def print_result(algo_name, clf, expr_te, label_te):
	label_predicted = clf.predict(expr_te)
	[sens, spec, accu] = calculate_confusion_matrix(label_te, label_predicted)
	sys.stdout.write('%s\t%.3f\t%.3f\t%.3f\t%.3f\n' % 
					(algo_name, sens, spec, accu, (sens+spec)/2))


def print_result2(algo_name, clf, expr_te, label_te):
	label_predicted = clf.predict(expr_te)
	# [sens, spec] = calculate_poly_detection(label_te, label_predicted)
	[sens, spec] = calculate_cancer_detection(label_te, label_predicted)
	sys.stdout.write('%s\t%.3f\t%.3f\t%.3f\t%.3f\n' % 
					(algo_name, sens, spec, 0, (sens+spec)/2))


## Arguments
N = int(sys.argv[1])
S = str(sys.argv[2]) ## True to load the N evaluation sets
S = True if S.lower() == 'true' else False
if S: 
	sys.stderr.write("Loading existed json\n")
else:
	sys.stderr.write("Generating new json\n")

## Project directories
dir_proj = '/Users/KANG/geneoscopy_dev/data/run_proj_batch1-17_3/'
# file_expr = dir_proj + '/training/training_set.txt'
file_expr = dir_proj + '/training/training_set_full.txt'

T = 200 #num of top genes

#### Evaluation of machine learning models ####
cross_valid = False

## Load training set
expr_full = np.loadtxt(file_expr, dtype=str, delimiter='\t')
feat_header = expr_full[0,:]
samples = expr_full[1:,0]
labels = expr_full[1:,1]

## Randomly select N times
label_keys = ['C', 'P', 'N']

if S:
	with open(dir_proj + '/validation/random_split_indices.json') as json_reader:
		label_dict = json.load(json_reader)
	for k in label_dict.keys(): ## replace the str key with numberic key
		label_dict[int(k)] = label_dict[k]
		del label_dict[k]
else:
	label_dict = {}


os.system('rm -rf tmp/; mkdir tmp/')

for i in range(N):
	sys.stdout.write('##### Set '+ str(i+1) + ' #####\n')
	sys.stderr.write('##### Set '+ str(i+1) + ' #####\n')

	if not S:
		## Randomly split
		label_dict[i] = {'training':{}, 'testing':{}}
		for k in label_keys:
			indx = np.where(labels == k)[0]
			indx_tr = np.random.choice(indx, size=np.floor(len(indx)*.8), replace=False)
			indx_te = np.setdiff1d(indx, indx_tr)
			label_dict[i]['training'][k] = list(indx_tr)
			label_dict[i]['testing'][k] = list(indx_te)

	## LIMMA wrapper
	# write_valid_chips(label_dict[i]['training'], samples, labels, 'tmp/valid_chips_'+str(i)+'.txt')

	# while not os.path.isfile('tmp/valid_chips_'+str(i)+'.txt'):
	# 	time.sleep(1)
	# os.system('Rscript ../../../scripts/de_analysis.r ../training/chipdata.txt tmp/valid_chips_'+str(i)+'.txt N_vs_P_vs_C 1 0 tmp/top_de_genes_'+str(i)+'.txt; cp tmp/top_de_genes_'+str(i)+'.txt tmp/tmp.txt; head -'+ str(T+1) +' tmp/tmp.txt > tmp/top_de_genes_'+str(i)+'.txt; rm tmp/tmp.txt')

	# while not os.path.isfile('tmp/top_de_genes_'+str(i)+'.txt'):
	# 	time.sleep(1)
	# de_genes = load_de_genes('tmp/top_de_genes_'+str(i)+'.txt')
	de_genes = load_de_genes('../training/top_de_genes.txt')
	expr = filter_features(expr_full, de_genes)


	## Random forest 

	if cross_valid:
		## sklearn model selection
		from sklearn.model_selection import GridSearchCV
		rf = RandomForestClassifier()
		hyperparams = {'n_estimators': [250, 500, 1000],
						'criterion': ['gini', 'entropy'],
						'class_weight': [None, 'balanced']}
		clf = GridSearchCV(rf, hyperparams, cv=parsed.cross_valid, n_jobs=4)
		clf.fit(expr_tr, label_tr)
		params = parse_cv_result(clf)
	else:
		params = {'n_estimators': 25,
						'criterion': 'gini',
						'class_weight': None}
	

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = RandomForestClassifier(n_estimators=params['n_estimators'], 
									criterion=params['criterion'],
									class_weight=params['class_weight'],
									oob_score=True,
									n_jobs=4, 
									verbose=False)
	clf.fit(expr_tr, label_tr)

	## test the model
	# print_result("Random forest", clf, expr_te, label_te)
	print_result2("Random forest", clf, expr_te, label_te)


	## C-SVM

	# if cross_valid:
	# 	## sklearn model selection
	# 	from sklearn.model_selection import GridSearchCV
	# 	svm = SVC()
		
	# 	# hyperparams = {'C': [.5, 1., 1.5, 2., 3,4,5,8,10],
	# 	# 				'kernel':['rbf'],
	# 	# 				# 'kernel': ['rbf', 'linear', 'poly', 'sigmoid'],
	# 	# 				'class_weight': [None]}
	# 	# clf = GridSearchCV(svm, hyperparams, cv=parsed.cross_valid, n_jobs=4)

	# 	from sklearn.model_selection import RandomizedSearchCV
	# 	import scipy.stats as ss
	# 	hyperparams = {'C': ss.expon(scale=10), #randomized parameters
	# 					'kernel':['rbf'],
	# 					# 'kernel': ['rbf', 'linear', 'poly', 'sigmoid'],
	# 					'class_weight': [None]}
	# 	clf = RandomizedSearchCV(svm, hyperparams, n_iter=500, cv=parsed.cross_valid, n_jobs=4)
		
	# 	clf.fit(expr_tr, label_tr)
	# 	params = parse_cv_result(clf)
	# else:
	# 	params = {'C': 1.1, #1.1 for 330 samples, 4.5 for 273 samples
	# 				'kernel': 'rbf',
	# 				'class_weight': None}

	# ## prepare data
	# label_split = label_dict[i]
	# [expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	# ## train the model
	# clf = SVC(C=params['C'], 
	# 			kernel=params['kernel'], 
	# 			class_weight=params['class_weight'],
	# 			probability=True, 
	# 			verbose=False)
	# clf.fit(expr_tr, label_tr)

	# ## test the model
	# print_result("C-SVM", clf, expr_te, label_te)
	# print_result2("C-SVM", clf, expr_te, label_te)


	## Nu-SVM

	if cross_valid:
		## sklearn model selection
		from sklearn.model_selection import GridSearchCV
		svm = NuSVC()
		
		# hyperparams = {'nu': [.5, 1., 1.5, 2., 3,4,5,8,10],
		# 				'kernel':['rbf'],
		# 				# 'kernel': ['rbf', 'linear', 'poly', 'sigmoid'],
		# 				'class_weight': [None]}
		# clf = GridSearchCV(svm, hyperparams, cv=parsed.cross_valid, n_jobs=4)

		from sklearn.model_selection import RandomizedSearchCV
		import scipy.stats as ss
		hyperparams = {'nu': ss.expon(scale=10), #randomized parameters
						'kernel':['rbf'],
						# 'kernel': ['rbf', 'linear', 'poly', 'sigmoid'],
						'class_weight': [None]}
		clf = RandomizedSearchCV(nusvm, hyperparams, n_iter=500, cv=parsed.cross_valid, n_jobs=4)
		
		clf.fit(expr_tr, label_tr)
		params = parse_cv_result(clf)
	else:
		params = {'nu': .82,
					'kernel': 'rbf',
					'class_weight': 'balanced'}

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = NuSVC(nu=params['nu'], 
				kernel=params['kernel'], 
				class_weight=params['class_weight'],
				probability=True, 
				verbose=False)
	clf.fit(expr_tr, label_tr)

	## test the model
	# print_result("Nu-SVM", clf, expr_te, label_te)
	print_result2("Nu-SVM", clf, expr_te, label_te)




	## Gradient boosting

	if cross_valid:
		## sklearn model selection
		from sklearn.model_selection import GridSearchCV
		gb = GradientBoostingClassifier()
		hyperparams = {'learning_rate': [.01, .0075, .005, .001, .0005], 
						'max_depth': [3],
						'subsample': [1, .8, .5],
						'n_estimators': [1000]}
		clf = GridSearchCV(gb, hyperparams, cv=parsed.cross_valid, n_jobs=4)
		clf.fit(expr_tr, label_tr)
		params = parse_cv_result(clf)
	else:
		params = {'learning_rate': 1, 
					'max_depth': 5,
					'subsample': 1.0,
					'n_estimators': 50}

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = GradientBoostingClassifier(learning_rate=params['learning_rate'], 
										n_estimators=params['n_estimators'], 
										max_depth=params['max_depth'], 
										subsample=params['subsample'], 
										verbose=False)
	clf.fit(expr_tr, label_tr)

	#test the model
	# print_result("Gradient boost", clf, expr_te, label_te)
	print_result2("Gradient boost", clf, expr_te, label_te)




	## AdaBoost

	if cross_valid:
		## sklearn model selection
		from sklearn.model_selection import GridSearchCV
		ab = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3))
		hyperparams = {'learning_rate': [.01, .0075, .005, .001, .0005], 
						'n_estimators': [1000]}
		clf = GridSearchCV(ab, hyperparams, cv=parsed.cross_valid, n_jobs=4)
		clf.fit(expr_tr, label_tr)
		params = parse_cv_result(clf)
	else:
		params = {'learning_rate': 1, 
					'n_estimators': 50}

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), 
								learning_rate=params['learning_rate'],
								n_estimators=params['n_estimators'])
	clf.fit(expr_tr, label_tr)

	#test the model
	# print_result("AdaBoost", clf, expr_te, label_te)
	print_result2("AdaBoost", clf, expr_te, label_te)




	# ## Neural network

	# ## prepare data
	# label_split = label_dict[i]
	# [expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	# ## train the model
	# clf = MLPClassifier(solver='adam', alpha=.0001, 
	# 					hidden_layer_sizes=(5, 2))
	# clf.fit(expr_tr, label_tr)

	# #test the model
	# print_result("Neural net", clf, expr_te, label_te)
	# # print_result2("Neural net", clf, expr_te, label_te)




	# ## Gaussian process

	# if cross_valid:
	# 	## sklearn model selection
	# 	from sklearn.model_selection import GridSearchCV
	# 	gb = GaussianProcessClassifier()
	# 	hyperparams = {}
	# 	clf = GridSearchCV(gb, hyperparams, cv=parsed.cross_valid, n_jobs=4)
	# 	clf.fit(expr_tr, label_tr)
	# 	params = parse_cv_result(clf)
	# else:
	# 	params = {}

	# ## prepare data
	# label_split = label_dict[i]
	# [expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	# ## train the model
	# clf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=1.0), 
	# 								optimizer="fmin_l_bfgs_b")
	# clf.fit(expr_tr, label_tr)

	# #test the model
	# print_result("Gauss process", clf, expr_te, label_te)
	# print_result2("Gauss process", clf, expr_te, label_te)



	## Naive Bayes

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = GaussianNB()
	clf.fit(expr_tr, label_tr)

	#test the model
	# print_result("Naive Bayes", clf, expr_te, label_te)
	print_result2("Naive Bayes", clf, expr_te, label_te)




	## KNN

	if cross_valid:
		## sklearn model selection
		pass
	else:
		params = {}

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = KNeighborsClassifier(n_neighbors=3, algorithm='ball_tree')
	clf.fit(expr_tr, label_tr)

	#test the model
	# print_result("kNN", clf, expr_te, label_te)
	print_result2("kNN", clf, expr_te, label_te)



	## Decision Tree

	if cross_valid:
		## sklearn model selection
		pass
	else:
		params = {}

	## prepare data
	label_split = label_dict[i]
	[expr_tr, label_tr, expr_te, label_te] = split_expr(expr, label_split)
	
	## train the model
	clf = DecisionTreeClassifier(max_depth=None)
	clf.fit(expr_tr, label_tr)

	#test the model
	# print_result("Decision tree", clf, expr_te, label_te)
	print_result2("Decision tree", clf, expr_te, label_te)



## Dump json data
# if not S:
# 	with open(dir_proj + '/validation/random_split_indices.json', 'w') as writer:
# 		json.dump(label_dict, writer)


