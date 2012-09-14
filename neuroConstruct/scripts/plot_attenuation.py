#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates and compares attenuation for the detailed and reduced
morphology models of the Golgi cell.

Attenuation is defined in the following way: stimulate a dendrite with
an aEPSP at a given distance from the soma, and measure the voltage
response at the injection site and at the soma. The attenuation is the
ratio between the peak of the somatic and the dendritic response.
"""
import sys
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
detailed_seg_ids = [sys.argv[2], sys.argv[4], sys.argv[6]]
detailed_locs = np.array([sys.argv[3], sys.argv[5], sys.argv[7]])
reduced_seg_ids = [4,5,6]
reduced_locs = np.array([35., 105., 175.])

data = {'reduced':[], 'detailed':[]}

fig, ax = plt.subplots()
colors = 'rgb'

for distance_index in [0, 1, 2]:
    seg_id_reduced = reduced_seg_ids[distance_index]
    seg_id_detailed = detailed_seg_ids[distance_index]
    color = colors[distance_index]
    sim_ref = timestamp + '_' + str(distance_index)
    sim_path = '../simulations/' + sim_ref
    data['reduced'].append([reduced_locs[distance_index],
			    np.loadtxt('{0}/Golgi_reduced_0.dat'.format(sim_path)),
			    np.loadtxt('{0}/Golgi_reduced_0.{1}.dat'.format(sim_path,
									    seg_id_reduced))])
    data['detailed'].append([detailed_locs[distance_index],
			     np.loadtxt('{0}/Golgi_Vervaeke_0.dat'.format(sim_path)),
			     np.loadtxt('{0}/Golgi_Vervaeke_0.{1}.dat'.format(sim_path,
									      seg_id_detailed))])

maxima = {'reduced':{}, 'detailed':{}}
baselines = {'reduced':{}, 'detailed':{}}
amplitudes = {'reduced':{}, 'detailed':{}}
attenuation = {}

for cell_type in data.keys():
    maxima[cell_type]['soma'] = np.array([data[cell_type][dist_index][1].max() for dist_index in [0,1,2]])
    maxima[cell_type]['dend'] = np.array([data[cell_type][dist_index][2].max() for dist_index in [0,1,2]])
    baselines[cell_type]['soma'] = np.array([data[cell_type][dist_index][1][-1] for dist_index in [0,1,2]])
    baselines[cell_type]['dend'] = np.array([data[cell_type][dist_index][2][-1] for dist_index in [0,1,2]])
    amplitudes[cell_type]['soma'] = maxima[cell_type]['soma'] - baselines[cell_type]['soma']
    amplitudes[cell_type]['dend'] = maxima[cell_type]['dend'] - baselines[cell_type]['dend']

    attenuation[cell_type] = amplitudes[cell_type]['soma']/amplitudes[cell_type]['dend']
    ax.plot([data[cell_type][dist_index][0] for dist_index in [0,1,2]], attenuation[cell_type],
	    marker='o', label=cell_type)

ax.legend(loc='best')
print attenuation



plt.show()
