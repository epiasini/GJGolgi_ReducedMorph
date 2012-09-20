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
    threshold = 'min20'
    for amplitude in current_amplitude_range:
	filename = '../simulations/{0}_{1}/Golgi_{2}_0.SPIKES_{3}.spike'.format(timestamp,
										amplitude,
										cell_name,
										threshold)
	try:
	    data.append(len(np.loadtxt(filename)))
	except TypeError:
	    data.append(0)
    ax.plot(current_amplitude_range, data, label=cell_name, marker='o')
ax.legend(loc='best')
plt.show()
