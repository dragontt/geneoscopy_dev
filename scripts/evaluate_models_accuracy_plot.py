#!/usr/bin/python
import numpy as np
import sys
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def customize_boxplot(plt, data, pos, facecolor):
	bp = plt.boxplot(data, positions=pos, whis=1, patch_artist=True)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='black', marker='*')

	colors = [facecolor] * len(data)
	for patch, color in zip(bp['boxes'], colors):
	    patch.set_facecolor(color)


filename = sys.argv[1]
T = int(sys.argv[2])


flag_mean_sens_spec = False


algos = ['Random forest', 'SVM', 'Gradient boosting', 'AdaBoost', 'Naive Bayes', 'kNN', 'Decision tree']
if flag_mean_sens_spec:
	colors = ["#FF6326", "#D45D3F", "#AA5859", "#7F5372", "#554D8C", "#2A48A5", "#0043BF"]
else:
	colors = ["#FF6326", "#D45D3F", "#AA5859", "#7F5372", "#554D8C", "#2A48A5", "#0043BF", "#00379F", "#002C7F", "#00215F", "#00163F", "#000B1F", "#000000"]
algo_indx = [1,2,0,3,4,6,5]

# algos = ['Random forest', 'SVM', 'AdaBoost', 'Gaussian process', 'kNN']
# # colors = ["#FF6326", "#BF5B4C", "#7F5372", "#3F4B98", "#0043BF"]
# colors = ["#FF6326", "#DF5F39", "#BF5B4C", "#9F575F", "#7F5372", "#5F4F85", "#3F4B98", "#1F47AB", "#0043BF"]
# algo_indx = [1,0,2,4,3]

N = len(algos)
# algo_indx = range(N)

## count frequencies
sens = np.zeros((N,T)) 
spec = np.zeros((N,T)) 
accu = np.zeros((N,T))
data = np.loadtxt(filename, delimiter='\t', usecols=[1,2,3,4])
for t in range(T):
	scores = data[t*N:(t+1)*N]
	sens[:,t] = scores[:,0]
	spec[:,t] = scores[:,1]
	accu[:,t] = scores[:,2]

## plot stack frequencies
indx = [i*3 for i in range(N)]
fig, ax = plt.subplots(figsize=(6, 5), dpi=150)

customize_boxplot(plt, sens[algo_indx,:].T, np.arange(0,(N*3), 3)+.5, '#507cb2')
customize_boxplot(plt, spec[algo_indx,:].T, np.arange(0,(N*3), 3)+1.5, '#73aa53')

plt.xticks(indx, np.array(algos)[algo_indx], rotation=30)
plt.ylim([0, 1])
plt.gcf().subplots_adjust(bottom=0.25)
# plt.show()

plt.savefig('results_'+ str(T) +'_runs.sens_spec.jpg', fmt='jpg')
# plt.savefig('results_'+ str(T) +'_runs.sens_spec.cancer_detection.jpg', fmt='jpg')
plt.close()




