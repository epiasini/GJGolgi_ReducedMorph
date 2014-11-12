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
figsize=(3.5,3)

timestamp = sys.argv[1]
current_amplitude_range = np.array([-25.e3, -20.e3, -15.e3, -10.e3, -5.e3, 0., 50.e3, 100.e3, 150.e3, 200.e3, 250.e3]) # (fA)     #np.arange(-25000, 500, 500)
sim_duration = 6000. #(ms)
transient = 1000. #(ms)

colors = ['k', 'r', 'g']

fig, ax = plt.subplots(figsize=figsize)

for k, cell_name in enumerate(['Solinas', 'Vervaeke', 'reduced']):
    data = []
    threshold = 'min20'
    for amplitude in current_amplitude_range:
	filename = '../simulations/{0}_{1:.0f}/Golgi_{2}_0.SPIKES_{3}.spike'.format(timestamp,
										amplitude,
										cell_name,
										threshold)
	try:
            spikes = np.loadtxt(filename)
            rate = 1000 * np.sum(spikes > transient)/(sim_duration - transient) # (Hz)
	    data.append(rate)
	except TypeError:
	    data.append(0)
    ax.plot(current_amplitude_range/1000.,
            data,
            label=cell_name,
            marker='o',
            markersize=5,
            #color=colors[k],
            linestyle='-')
ax.legend(loc='lower right')
ax.set_xlabel('Injected current (pA)')
ax.set_ylabel('Firing rate (Hz)')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.locator_params(tight=True, nbins=5)
plt.tight_layout()
fig.savefig('fig/fi.pdf')
plt.show()
