#! /usr/bin/env python
# -*- coding: utf-8 -*-
## plotting script for F-I analysis of single cell.
## Usage: pythonv fi_plot.py b1415721979.9
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import seaborn as sns

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)
figsize = (1.75,1.25) # (3.5,3) for full-size figure

timestamp = sys.argv[1]
current_amplitude_range =  np.arange(-25000, 500, 500) # (fA)   # np.array([-25.e3, -20.e3, -15.e3, -10.e3, -5.e3, 0., 50.e3, 100.e3, 150.e3, 200.e3, 250.e3])
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
            linewidth=1.5, # 2
            marker='', # 'o' for full-size figure
            markersize=3, # 5 for full-size figure
            #color=colors[k],
            linestyle='-')
#ax.legend(loc='lower right')
ax.set_xlabel('Injected current (pA)')
ax.set_ylabel('Firing rate (Hz)')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.locator_params(tight=False, nbins=4) # nbins=5 probably better for full-size figure
ax.set_xlim((-27, 2)) #ax.set_xlim((-35, 260))
ax.set_ylim((-1, 10)) #ax.set_ylim((-2, 60))
plt.tight_layout()
fig.savefig('fig/fi_inset.pdf')
plt.show()
