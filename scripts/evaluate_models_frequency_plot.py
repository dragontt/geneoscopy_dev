#!/usr/bin/python
import numpy as np
import sys
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


filename = sys.argv[1]
T = int(sys.argv[2])


flag_mean_sens_spec = False


# algos = ['Random forest', 'SVM', 'Gradient boosting', 'AdaBoost', 'Naive Bayes', 'kNN', 'Decision tree']
# if flag_mean_sens_spec:
# 	colors = ["#FF6326", "#D45D3F", "#AA5859", "#7F5372", "#554D8C", "#2A48A5", "#0043BF"]
# else:
# 	colors = ["#FF6326", "#D45D3F", "#AA5859", "#7F5372", "#554D8C", "#2A48A5", "#0043BF", "#00379F", "#002C7F", "#00215F", "#00163F", "#000B1F", "#000000"]
# algo_indx = [1,2,0,3,4,6,5]


algos = [.76, .78, .8, .82, .84, .86, .88, .9, .92]
algos = [str(a) for a in algos]
colors = ["#FF6326", "#DF5F39", "#BF5B4C", "#9F575F", "#7F5372", "#5F4F85", "#3F4B98", "#1F47AB", "#0043BF", "#003AA7", "#00328F", "#002977", "#00215F", "#001947", "#00102F", "#000817", "#000000"]
algo_indx = range(len(algos))


N = len(algos)


if flag_mean_sens_spec:
	## count frequencies
	freqs = np.zeros((N,N), dtype='int') 

	data = np.loadtxt(filename, delimiter='\t', usecols=[1,2,3,4])
	for t in range(T):
		scores = data[t*N:(t+1)*N]
		rnks = np.zeros(scores.shape, dtype='int')
		for j in range(4):
			rnks[:,j] = np.ones(N)*N - rankdata(scores[:,j], method='max') 

		for n in range(N):
			freqs[n,rnks[n,3]] += 1

	freqs = np.array(freqs)/float(T)

	## plot stack frequencies
	indx = .1+np.arange(N)
	bottom = np.zeros(N)

	plt.figure(num=None, figsize=(6, 6), dpi=150)
	for i in reversed(range(N)):
		plt.bar(indx, freqs[algo_indx,i], .9, color=colors[i], bottom=bottom)
		bottom += freqs[algo_indx,i]
	plt.xticks(indx, np.array(algos)[algo_indx], rotation=90)
	plt.gcf().subplots_adjust(bottom=0.25)
	# plt.show()

else:
	## count frequencies
	freqs = np.zeros((N,(N*2-1)), dtype='int') 

	data = np.loadtxt(filename, delimiter='\t', usecols=[1,2,3,4])
	for t in range(T):
		scores = data[t*N:(t+1)*N]
		rnks = np.zeros(scores.shape, dtype='int')
		for j in range(4):
			rnks[:,j] = np.ones(N)*N - rankdata(scores[:,j], method='max') 

		for n in range(N):
			x = rnks[n,0] + rnks[n,1]
			freqs[n,x] += 1

	freqs = np.array(freqs)/float(T)

	## plot stack frequencies
	indx = .1+np.arange(N)
	bottom = np.zeros(N)

	plt.figure(num=None, figsize=(6, 4), dpi=300)
	for i in reversed(range(N*2-1)):
		plt.bar(indx, freqs[algo_indx,i], .9, color=colors[i], bottom=bottom)
		bottom += freqs[algo_indx,i]
	algos_print = ['\n'.join(a.split(' ')) for a in algos]
	plt.xticks(indx+.5, np.array(algos_print)[algo_indx], rotation=40)
	plt.tick_params(axis=u'x', which=u'both',length=0)
	plt.ylim([0,1])
	plt.gcf().subplots_adjust(bottom=0.25)
	# plt.show()
	# """

# plt.savefig('results_'+ str(T) +'_runs.jpg', fmt='jpg')
# plt.savefig('results_'+ str(T) +'_runs.cancer_detection.jpg', fmt='jpg')
plt.savefig('results_'+ str(T) +'_runs.nusvm_hyperparams.jpg', fmt='jpg')
plt.close()

# """
## make colorbar
plt.figure(num=None, figsize=(.5, 4), dpi=300)
fig, ax = plt.subplots(figsize=(.5, 4), dpi=300)
for i in range(N*2-1):
	ax.bar(0, 10./9, 1, color=colors[(N-1)*2-i], bottom=10./9*i)
plt.tick_params(axis=u'both', which=u'both',length=0)
plt.ylim([0,14.45])
ax.set_xticklabels('')
ax.yaxis.tick_right()
ax.yaxis.set_major_formatter(ticker.NullFormatter())
ax.yaxis.set_minor_locator(ticker.FixedLocator([i*2.25+.5 for i in range(N*2-1)]))
ax.yaxis.set_minor_formatter(ticker.FixedFormatter(range(1,8)[::-1]))
plt.gcf().subplots_adjust(right=0.6, bottom=0.25)
# plt.show()
# plt.savefig('results_'+ str(T) +'_runs.colorbar.jpg', fmt='jpg')
plt.close()
# """


