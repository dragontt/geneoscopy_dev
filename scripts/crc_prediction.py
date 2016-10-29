#/usr/bin/python
import sys
import argparse
import numpy
from sklearn.externals import joblib

learning_algorithms = ['random_forest', 'svm', 'neural_net', 'grad_boosting']

def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-a', '--learning_algorithm', dest='learning_algorithm', default='random_forest', help='options: %s' % learning_algorithms)
	parser.add_argument('-i', '--input_filename', dest='input_filename')
	parser.add_argument('-m', '--model_filename', dest='model_filename')
	parsed = parser.parse_args(argv[1:])
	return parsed

def parse_data(filename, label_col, data_col_start, data_col_end=10000):
	data = numpy.loadtxt(filename, dtype=str, delimiter='\t')
	gene_id = data[0, data_col_start:data_col_end]
	sample_id = data[1:, 0]
	expr = numpy.array(data[1:, data_col_start:data_col_end], dtype=numpy.float32)
	label = data[1:, label_col]
	return [gene_id, sample_id, expr, label]

def main(argv):
	# parse data
	parsed = parse_args(argv)
	[gene_id, sample_id, expr_te, label_te] = parse_data(parsed.input_filename, 1, 2)

	label_unique= numpy.unique(label_te)
	label_count = numpy.array([len(numpy.where(label_te == l)[0]) for l in label_unique])

	print "Validation set dimension:", expr_te.shape[0], "samples x", expr_te.shape[1], "features"
	print "CRC labels:", label_unique, ", counts:", label_count 


	##### Random Forest #####
	if parsed.learning_algorithm.lower() ==	'random_forest':
		from sklearn.ensemble import RandomForestClassifier
		
		# predict on validation set
		clf = joblib.load(parsed.model_filename)
		label_predicted = clf.predict(expr_te)
		probability_predicted = clf.predict_proba(expr_te)
		accuracy_predicted = clf.score(expr_te, label_te)
		summary = numpy.hstack((sample_id[numpy.newaxis].T, label_te[numpy.newaxis].T, label_predicted[numpy.newaxis].T, probability_predicted))

		# print messages
		print "sample_id true_label predict_label", clf.classes_
		print '\n'.join(' '.join(str(cell) for cell in row) for row in summary)
		print "Prediction accuracy:", accuracy_predicted


	##### SVM #####
	elif parsed.learning_algorithm.lower() == 'svm':
		from sklearn.svm import SVC

		# predict on validation set
		clf = joblib.load(parsed.model_filename)
		label_predicted = clf.predict(expr_te)
		probability_predicted = clf.predict_proba(expr_te)
		accuracy_predicted = clf.score(expr_te, label_te)
		summary = numpy.hstack((sample_id[numpy.newaxis].T, label_te[numpy.newaxis].T, label_predicted[numpy.newaxis].T, probability_predicted))

		# print message 
		print "sample_id true_label predict_label", clf.classes_
		print '\n'.join(' '.join(str(cell) for cell in row) for row in summary)
		print "Prediction accuracy:", accuracy_predicted


	##### Neural Network #####
	elif parsed.learning_algorithm.lower() == 'neural_net':
		from sklearn.neural_network import BernoulliRBM

		# predict on validation set
		clf = joblib.load(parsed.model_filename)
		label_predicted = clf.predict(expr_te)
		accuracy_predicted = len([label_predicted[i] for i in range(len(label_predicted)) if (label_predicted[i] == label_te[i])]) / float(len(label_predicted))

		# print message 
		print "Prediction accuracy:", accuracy_predicted


	##### Gradient Boosting #####
	elif parsed.learning_algorithm.lower() == 'grad_boosting':
		from sklearn.ensemble import GradientBoostingClassifier
		
		# predict on validation set
		clf = joblib.load(parsed.model_filename)
		label_predicted = clf.predict(expr_te)
		probability_predicted = clf.predict_proba(expr_te)
		accuracy_predicted = clf.score(expr_te, label_te)
		summary = numpy.hstack((sample_id[numpy.newaxis].T, label_te[numpy.newaxis].T, label_predicted[numpy.newaxis].T, probability_predicted))

		# print messages
		print "sample_id true_label predict_label", clf.classes_
		print '\n'.join(' '.join(str(cell) for cell in row) for row in summary)
		print "Prediction accuracy:", accuracy_predicted


	##### Gaussian Process #####
	elif parsed.learning_algorithm.lower() == 'gauss_process':
		from sklearn.gaussian_process import GaussianProcessClassifier
		from sklearn.gaussian_process.kernels import RBF
		
		# predict on validation set
		clf = joblib.load(parsed.model_filename)
		label_predicted = clf.predict(expr_te)
		probability_predicted = clf.predict_proba(expr_te)
		accuracy_predicted = clf.score(expr_te, label_te)
		summary = numpy.hstack((sample_id[numpy.newaxis].T, label_te[numpy.newaxis].T, label_predicted[numpy.newaxis].T, probability_predicted))

		# print messages
		print "sample_id true_label predict_label", clf.classes_
		print '\n'.join(' '.join(str(cell) for cell in row) for row in summary)
		print "Prediction accuracy:", accuracy_predicted

	else:
		sys.exit('Improper learning algorithm option given.')

if __name__ == "__main__":
    main(sys.argv)

