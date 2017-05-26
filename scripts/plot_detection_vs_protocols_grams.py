#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.axes.Axes import violinplot
from collections import OrderedDict


def parse_dicts(file):
	f = open(file, 'r')
	lines = f.readlines()
	f.close()
	labels = {}
	grams = {}
	protocols = {'1':[], '2':[], '3':[], '4':[]}
	for i in range(len(lines)):
		if lines[i].startswith('2'):
			line = lines[i].strip().split('\t')
			labels[line[0]] = line[3]
			grams[line[0]] = float(line[1]) if line[1] != '' else -1
			protocols[line[2]].append(line[0])
	return (labels, grams, protocols)


def parse_model_output(file):
	f = open(file, 'r')
	lines = f.readlines()
	f.close()
	detections = {'tp':[], 'fp':[], 'fn':[], 'tn':[]}
	for i in range(1,len(lines)):
		line = lines[i].strip().split('\t')
		if line[3] == '1':
			detections['tp'].append(line[0])
		elif line[4] == '1':
			detections['fp'].append(line[0])
		elif line[5] == '1':
			detections['fn'].append(line[0])
		else:
			detections['tn'].append(line[0])
	return detections


def plot_grams_vs_protocols(labels, grams, protocols, detections, tr_samples, file_fig):
	shifts = {'Cancer':0, 'Polyp':1, 'Normal':2}
	fig, ax = plt.subplots(figsize=(7,6), dpi=150)
	i = 0
	for p in sorted(protocols.keys()):
		for s in protocols[p]:		
			if s in tr_samples:
				y = grams[s]
				x = np.random.normal(i+shifts[labels[s]], .1, size=1)
				# ax.scatter(x, y, s=50, c='k', alpha=.3, marker='+', label='Training')
			else:
				if s in detections['tp']:
					y = grams[s]
					x = np.random.normal(i+shifts[labels[s]], .1, size=1)
					ax.scatter(x, y, s=50, c='b', alpha=.3, marker='o', label='True positive')
				elif s in detections['tn']:
					y = grams[s]
					x = np.random.normal(i+shifts[labels[s]], .1, size=1)
					ax.scatter(x, y, s=50, c='g', alpha=.3, marker='o', label='True negeative')
				elif s in detections['fp']:
					y = grams[s]
					x = np.random.normal(i+shifts[labels[s]], .1, size=1)
					ax.scatter(x, y, s=80, c='r', alpha=.8, marker='*', label='False positive')
				elif s in detections['fn']:
					y = grams[s]
					x = np.random.normal(i+shifts[labels[s]], .1, size=1)
					ax.scatter(x, y, s=80, c='y', alpha=.8, marker='*', label='False negative')
		i += 3
	plt.xlim([-.5,11.5])
	plt.ylim([0,25])
	ax.set_xlabel('\n\nExtraction protocol number/type')
        ax.set_ylabel('Stool grams')
	ax.set_xticks(range(12))
	ax.set_xticklabels(['Cancer','Adenoma','Normal']*4, rotation=90)
	h, l = plt.gca().get_legend_handles_labels()
	bl = OrderedDict(zip(l, h))
	bl_values = bl.values()
	bl_keys = bl.keys()
	bl_values[0], bl_values[1] = bl.values()[1], bl.values()[0]
	bl_keys[0], bl_keys[1] = bl.keys()[1], bl.keys()[0]
	plt.legend(bl_values[::-1], bl_keys[::-1])
	plt.tight_layout()
	plt.savefig(file_fig, fmt='svg')


def plot_grams_vs_detections(labels, grams, protocols, detections, tr_samples, file_fig):
	fig, ax = plt.subplots(figsize=(6,4), dpi=150)
	data = []
	data.append([grams[s] for s in detections['tp']])
	data.append([grams[s] for s in detections['tn']])
	data.append([grams[s] for s in detections['fp']])
	data.append([grams[s] for s in detections['fn']])
	data.append([grams[s] for s in tr_samples])
	ax.violinplot(data, showmedians=True)
	ax.set_xticks(range(1,6))
	# ax.set_xticklabels(['True positive', 'True nagative', 'False positive', 'False negative', 'Training'], rotation=30)
	ax.set_xticklabels(['True positive', 'True nagative', 'False positive', 'False negative'], rotation=30)
	ax.set_ylabel('Gram')
	plt.tight_layout()
	plt.savefig(file_fig, fmt='svg')




labels, grams, protocols = parse_dicts('sample_inv_lite.txt')
detections = parse_model_output('results_04232017.txt')
tr_samples = np.loadtxt('../../training/valid_chips.txt', dtype=str, usecols=[0], delimiter='.')

plot_grams_vs_protocols(labels, grams, protocols, detections, tr_samples, 'plot_grams_vs_protocols.svg')
# plot_grams_vs_detections(labels, grams, protocols, detections, tr_samples, 'plot_grams_vs_detection.svg')
