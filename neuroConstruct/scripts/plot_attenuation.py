#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

target_segment_index_range = [3]
target_fraction_along_range = [.1,.2,.3,.4,.5,.6,.7,.8,.9]

timestamp = 1347010416.91

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
