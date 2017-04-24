#!/usr/bin/python
import numpy as np
import json
import os
import time
from scipy.stats import rankdata

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.svm import NuSVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF


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


## Project directories
dir_proj = '/Users/KANG/geneoscopy_dev/data/run_proj_batch1-17_3/'
# file_expr = dir_proj + '/training/training_set.txt'
file_expr = dir_proj + '/training/training_set_full.txt'
N = 10
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
label_dict = {}


hyperparams = [.5, .6, .7, .82, .9]
eval_metrics = np.zeros((10, len(hyperparams)))


for i in range(N):
	print '*** Set '+ str(i+1) + ' ***'

	## Randomly split
	label_dict[i] = {'training':{}, 'testing':{}}
	for k in label_keys:
		indx = np.where(labels == k)[0]
		indx_tr = np.random.choice(indx, size=np.floor(len(indx)*.8), replace=False)
		indx_te = np.setdiff1d(indx, indx_tr)
		label_dict[i]['training'][k] = list(indx_tr)
		label_dict[i]['testing'][k] = list(indx_te)

	## LIMMA wrapper
	write_valid_chips(label_dict[i]['training'], samples, labels, 'valid_chips_'+str(i)+'.txt')

	while not os.path.isfile('valid_chips_'+str(i)+'.txt'):
		time.sleep(1)
	print('analyzing DE genes for this set ...')
	os.system('Rscript ../../../scripts/de_analysis.r ../training/chipdata.txt valid_chips_'+str(i)+'.txt N_vs_P_vs_C 1 0 top_de_genes_'+str(i)+'.txt; cp top_de_genes_'+str(i)+'.txt tmp.txt; head -'+ str(T+1) +' tmp.txt > top_de_genes_'+str(i)+'.txt; rm tmp.txt')

	while not os.path.isfile('top_de_genes_'+str(i)+'.txt'):
		time.sleep(5)
	de_genes = load_de_genes('top_de_genes_'+str(i)+'.txt')
	expr = filter_features(expr_full, de_genes)


	## Nu-SVM
	for j in range(len(hyperparams)):
		p = hyperparams[j]

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
			params = {'nu': p,
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
		label_predicted = clf.predict(expr_te)
		[sens, spec, accu] = calculate_confusion_matrix(label_te, label_predicted)
		print 'Nu-SVM\t', p,'\t', sens, spec, accu


		eval_metrics[i, j] = accu

print ""
rankings = []
for i in range(eval_metrics.shape[0]):
	rankings.append( np.ones(len(hyperparams))*(len(hyperparams)+1) - rankdata(eval_metrics[i,:], method='max') )
rankings = np.array(rankings)
print rankings
print ""
print np.mean(rankings, axis=0)

## Dump json data
with open(dir_proj + '/validation/random_split_indices.json', 'w') as writer:
	json.dump(label_dict, writer)


