#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
current_amplitude_range = [-10, -5, 0, 5, 10, 15, 20]

fig, ax = plt.subplots()

for cell_name in ['reduced', 'Vervaeke', 'Solinas']:
    data = []
    for amplitude in current_amplitude_range:
	filename = '../simulations/{0}_{1}/Golgi_{2}_0.SPIKES_min20.spike'.format(timestamp,
										  amplitude,
										  cell_name)
	data.append(len(np.loadtxt(filename)))
    ax.plot(current_amplitude_range, data, label=cell_name, marker='o')
ax.legend(loc='best')
plt.show()
