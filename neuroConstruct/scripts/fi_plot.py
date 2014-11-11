#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import seaborn as sns

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{euler}']
matplotlib.rcParams['font.sans-serif'].insert(0, 'Bitstream Vera Sans')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.size'] = 8.0
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['xtick.labelsize'] = 'medium'
matplotlib.rcParams['ytick.labelsize'] = 'medium'
matplotlib.rcParams['axes.labelsize'] = 'medium'
figsize=(1.75,1.25)

timestamp = sys.argv[1]
current_amplitude_range = np.arange(-25000, 500, 500)
sim_duration = 6. #(s)
transient = 1. #(s)

colors = ['k', 'r', 'g']

fig, ax = plt.subplots(figsize=figsize)

for k, cell_name in enumerate(['Solinas', 'Vervaeke', 'reduced']):
    data = []
    threshold = 'min20'
    for amplitude in current_amplitude_range:
	filename = '../simulations/{0}_{1}/Golgi_{2}_0.SPIKES_{3}.spike'.format(timestamp,
										amplitude,
										cell_name,
										threshold)
	try:
            spikes = np.loadtxt(filename)
            rate = np.sum(spikes > transient)/(sim_duration - transient)
	    data.append(rate)
	except TypeError:
	    data.append(0)
    ax.plot(current_amplitude_range/1000.,
            data,
            label=cell_name,
            marker='o',
            markersize=3,
            #color=colors[k],
            linestyle='')
#ax.legend(loc='lower right')
ax.set_xlabel('Injected current (pA)')
ax.set_ylabel('Firing rate (Hz)')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.locator_params(tight=True, nbins=4)
plt.tight_layout()
fig.savefig('test/fi.pdf')
plt.show()
