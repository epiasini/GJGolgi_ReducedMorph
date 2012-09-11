#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

if len(sys.argv) > 1:
    timestamp = sys.argv[1]
else:
    timestamp = 1347282961.43

target_fraction_along_range = [.1]

data = []

fig, ax = plt.subplots()

for cell_name in ['reduced', 'Vervaeke']:
    if cell_name == 'reduced': target_segments = [4,5,6]
    else: target_segments = [1526, 1545, 1646]
    for target_segment_index in target_segments:
	for target_fraction_along in target_fraction_along_range:
	    filename = '../simulations/{0}_{1}_{2}/Golgi_{3}_0.dat'.format(timestamp,
									       target_segment_index, target_fraction_along, cell_name)
	    data.append([target_segment_index,
			 target_fraction_along,
			 np.loadtxt(filename)])
	    ax.plot(data[-1][2], label='{0}_{1}'.format(target_segment_index,
							target_fraction_along))

ax.legend('off')
plt.show()
