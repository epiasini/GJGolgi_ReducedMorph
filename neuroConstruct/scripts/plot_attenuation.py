#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

if len(sys.argv) > 1:
    timestamp = sys.argv[1]
else:
    timestamp = 1347282961.43

target_segment_index_range = [3]
target_fraction_along_range = [.1,.5,.9]

data = []

fig, ax = plt.subplots()

for target_segment_index in target_segment_index_range:
    for target_fraction_along in target_fraction_along_range:
	filename = '../simulations/{0}_{1}_{2}/Golgi_reduced_0.dat'.format(timestamp,
									   target_segment_index,
									   target_fraction_along)
	data.append([target_segment_index, target_fraction_along, np.loadtxt(filename)])
	ax.plot(data[-1][2], label='{0}_{1}'.format(target_segment_index,
						    target_fraction_along))

ax.legend()
plt.show()
