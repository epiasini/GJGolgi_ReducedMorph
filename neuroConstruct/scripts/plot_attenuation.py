#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
id0 = sys.argv[2]
id1 = sys.argv[3]
id2 = sys.argv[4]

target_fraction_along_range = [.1]

data = []

fig, ax = plt.subplots()

for cell_name in ['reduced', 'Vervaeke']:
    if cell_name == 'reduced': target_segments = [4,5,6]
    else: target_segments = [id0, id1, id2]
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
