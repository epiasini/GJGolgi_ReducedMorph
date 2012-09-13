#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
ids = [sys.argv[2], sys.argv[4], sys.argv[6]]
vervaeke_locs = np.array([sys.argv[3], sys.argv[5], sys.argv[7]])
reduced_locs = np.array([9.18, 33.88, 124.69])

target_fraction_along_range = [.5]

data = []

fig, ax = plt.subplots()

for cell_name in ['reduced', 'Vervaeke']:
    if cell_name == 'reduced': target_segments = [4,5,6]
    else: target_segments = ids
    for target_segment_index in target_segments:
	for target_fraction_along in target_fraction_along_range:
	    filename = '../simulations/{0}_{1}_{2}/Golgi_{3}_0.dat'.format(timestamp,
									   target_segment_index, target_fraction_along, cell_name)
	    data.append([target_segment_index,
			 target_fraction_along,
			 np.loadtxt(filename)])
	    ax.plot(data[-1][2][6012:], label='{0}_{1}'.format(target_segment_index,
							target_fraction_along))

traces = np.array([x[2] for x in data])
maxima = traces.max(axis=1)
baselines = traces[:,-1]
peak_amplitudes = maxima - baselines
#print(maxima)
#print(baselines)
#print peak_amplitudes

reduced_amplitudes = peak_amplitudes[0:3]
vervaeke_amplitudes = peak_amplitudes[3:]

ampl_fig, ampl_ax = plt.subplots()
ampl_ax.plot(reduced_locs, reduced_amplitudes, marker='o')
ampl_ax.plot(vervaeke_locs, vervaeke_amplitudes, marker='o')

plt.show()
