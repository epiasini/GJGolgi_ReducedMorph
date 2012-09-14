#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
detailed_seg_ids = [sys.argv[2], sys.argv[4], sys.argv[6]]
detailed_locs = np.array([sys.argv[3], sys.argv[5], sys.argv[7]])
reduced_seg_ids = [4,5,6]
reduced_locs = np.array([9.18, 33.88, 124.69])

data = {'reduced':[], 'detailed':[]}

fig, ax = plt.subplots()

colors = 'rgb'

for distance_index in [0, 1, 2]:
    color = colors[distance_index]
    sim_ref = timestamp + '_' + str(distance_index)
    sim_path = '../simulations/' + sim_ref
    data['reduced'].append([reduced_locs[distance_index],
			    np.loadtxt('{0}/Golgi_reduced_0.dat'.format(sim_path))])
    data['detailed'].append([detailed_locs[distance_index],
			     np.loadtxt('{0}/Golgi_Vervaeke_0.dat'.format(sim_path))])
    ax.plot(data['reduced'][distance_index][1], linestyle='-', color=color)
    ax.plot(data['detailed'][distance_index][1], linestyle=':', color=color)

print(data)

maxima = {}
baselines = {}
amplitudes = {}
for cell_type in data.keys():
    print data[cell_type][0][1]
    maxima[cell_type] = np.array([data[cell_type][dist_index][1].max() for dist_index in [0,1,2]])
    baselines[cell_type] = np.array([data[cell_type][dist_index][1][-1] for dist_index in [0,1,2]])
    amplitudes[cell_type] = maxima[cell_type] - baselines[cell_type]

print amplitudes

plt.show()
