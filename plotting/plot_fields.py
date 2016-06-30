#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 8

plt.rc('font', size=16)

vecs=np.loadtxt("../output/pdfout.txt")
positions=np.loadtxt("../output/positions.txt")
print vecs.shape
print positions

xvec=vecs[:,0]
pdf_vec=vecs[:,1]

fig,ax=plt.subplots()       

ax.plot(xvec,pdf_vec,"-")



# the histogram of the data
n, bins, patches = ax.hist(positions, 500, normed=1, facecolor='green', alpha=0.75)

# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#ax.plot(bins, 'r--', linewidth=1)



#ax.set_yscale('log')
#ax.set_ylim([1e-1,1])
#ax.set_ylim([0,1.0])
#ax.set_xlim([0,1.65e4])
ax.set_xlabel("X")
ax.set_ylabel("PDF")
#ax.set_title("Antineutrino With Decay, With Matter")


# Now add the legend with some customizations.
#legend = ax.legend(loc='center right', shadow=False)
#legend = ax.legend(bbox_to_anchor=(0.9, 0.6), bbox_transform=plt.gcf().transFigure, shadow=False)
"""
legend = ax.legend(loc='center right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(24)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
"""

ax.tick_params(axis='x', which='major', labelsize=16)
ax.tick_params(axis='y', which='major', labelsize=16)

plt.show()
#print sp.AtmosphericNeutrinoOscillationProbability(1,1,100*param.GeV,param.PI,param)


