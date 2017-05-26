#/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt

##Color choice:
#507cb2 <- blue
#73aa53 <- green
#7f4d91 <- purple


filename = sys.argv[1]
color = sys.argv[2]

f = open(filename, 'r')
lines = f.readlines()
f.close()

labels = []
values = []
for i in range(1, (len(lines)-1)):
	line = lines[i].strip().split('\t')
	labels.append(line[0])
	values.append(float(line[2].strip('%')))


indx = np.arange(len(labels)) + .1
fig, ax = plt.subplots(figsize=(4, 2.5), dpi=150)

plt.bar(indx, values, .5, color='#'+color)

plt.ylabel('Sensitivity (%)')
plt.xticks(indx+.25, labels, rotation=40)
plt.tick_params(axis=u'x', which=u'both',length=0)
plt.ylim([0, 100])
plt.gcf().subplots_adjust(bottom=0.35, left=.2)
# plt.show()

rects = ax.patches
for rect, value in zip(rects, values):
    height = rect.get_height()
    annot_text = ax.text(rect.get_x() + rect.get_width()/2, height - 12, ('%d%%' % value), 
    	ha='center', va='bottom', color='white')
    annot_text.set_fontsize(9)

plt.savefig(filename.strip('.txt')+'.pdf', fmt='pdf')
