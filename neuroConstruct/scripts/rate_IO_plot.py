#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Network I/O plotting script.
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
stim_rate_range = range(1, 440, 40)

fig, ax = plt.subplots()

cell_types = ['reduced', 'Vervaeke', 'Solinas']
rates = {}
threshold = 'min20'

for cell_type in cell_types:
    rates[cell_type] = []
    for amplitude in stim_rate_range:
	filename = '../simulations/{0}_{1}/Golgi_{2}_0.SPIKES_{3}.spike'.format(timestamp,
										amplitude,
										cell_type,
										threshold)
	rates[cell_type].append(len(np.loadtxt(filename)))
    ax.plot(stim_rate_range, rates[cell_type], marker='o', label=cell_type)

ax.legend(loc='best')
ax.set_xlabel('stimulation rate (Hz)')
ax.set_ylabel('firing rate (Hz)')
plt.show()
