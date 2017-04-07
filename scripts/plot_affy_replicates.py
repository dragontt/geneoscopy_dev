#!/usr/bin/python
import sys
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm

def plot_replicates(reps_dict, sample_ids, data, dir_figures):
	for k in reps_dict.keys():
		rep_indx = []
		for rep in reps_dict[k]:
			rep_indx.append(np.where(sample_ids == rep)[0][0])

		##fit a line, calculate concordance, and plot
		model = sm.OLS(data[:,rep_indx[0]], data[:,rep_indx[1]])
		results = model.fit()
		# print results.summary()
		
		plt.figure(num=None, figsize=(5, 4.5), dpi=300)
		matplotlib.rcParams.update({'font.size': 16, 'font.family':'Arial'})  
		plt.scatter(data[:,rep_indx[0]], data[:,rep_indx[1]], alpha=.3, lw=0)
		plt.plot(data[:,rep_indx[0]], results.params[0]*data[:,rep_indx[0]], 'r')
		
		plt.text(15, 1, '$adj. R^2 = %.3f$' % results.rsquared_adj)
		plt.xlabel('Replicate 1')
		plt.ylabel('Replicate 2')
		plt.gcf().subplots_adjust(bottom=0.15)
		plt.savefig(dir_figures+'/'+k+'.jpg', format='jpg')


##replicate pairs
tech_reps = {'tech.P63': ['27998', '28002'],
			 'tech.C20-3': ['27999', '28003'],
			 'tech.H99': ['28000', '28004'],
			 'tech.H121': ['28001', '28005']}
bio_reps = {'bio.P63.A': ['27998', '27586'],
			'bio.P63.B': ['28002', '27586'],
			'bio.C20-3.A': ['27999', '27782'],
			'bio.C20-3.B': ['28003', '27782'],
			'bio.H99.A': ['28000', '27442'],
			'bio.H99.B': ['28004', '27442'],
			'bio.H121.A': ['28001', '27732'],
			'bio.H121.B': ['28005', '27732']}
freezer_reps = {'freezer.C36-4': ['27907', '27589'],
				'freezer.C34-2': ['27908', '27599'],
				'freezer.H70': ['27909', '27387'],
				'freezer.H92': ['27910', '27435'],
				'freezer.P98': ['27911', '27543'],
				'freezer.P105': ['27912', '27544']}

dir_project = '/Users/KANG/geneoscopy_dev/data/20170406_combined_projects_360samples/'
filename_expr = dir_project + 'gc-correction.scale-intensities.rma-bg.quant-norm.pm-only.med-polish.summary.txt'

##load data
data = np.loadtxt(filename_expr, dtype=str, delimiter='\t')
samples = data[0,1:]
genes = data[1:,0]
data = np.array(data[1:,1:], dtype=float)

##parse sample IDs
sample_ids = np.empty((len(samples),1), dtype='|S20')
for i in range(len(samples)):
	sample_ids[i] = re.split('\s|\.', samples[i])[0]

##plot replicates
plot_replicates(tech_reps, sample_ids, data, dir_project+'/figures')
plot_replicates(bio_reps, sample_ids, data, dir_project+'/figures')
plot_replicates(freezer_reps, sample_ids, data, dir_project+'/figures')

