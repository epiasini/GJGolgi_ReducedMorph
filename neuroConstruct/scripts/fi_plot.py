#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
current_amplitude_range = np.arange(-25000, 500, 500)
sim_duration = 6. #(s)
transient = 1. #(s)

colors = ['k', 'r', 'g']

fig, ax = plt.subplots()

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
            color=colors[k],
            linestyle='')
ax.legend(loc='best')
ax.set_xlabel('injected current (pA)')
ax.set_ylabel('firing rate (Hz)')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
fig.savefig('test.png')
plt.show()
