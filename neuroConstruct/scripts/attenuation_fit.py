#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates and compares attenuation for the detailed and reduced
morphology models of the Golgi cell.

Attenuation is defined in the following way: stimulate a dendrite with
an aEPSP at a given distance from the soma, and measure the voltage
response at the injection site and at the soma. The attenuation is the
ratio between the peak of the somatic and the dendritic response.

Note that this is _not_ how Vervaeke2012 defines dendritic attenuation.
"""
import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def attenuation_function(x, l):
    return np.exp(-x/l)

initial_a = 1.
initial_l = 46.

n_points = 18
timestamp = sys.argv[1]
detailed_seg_ids = [sys.argv[n] for n in range(2, 38, 2)]
detailed_locs = np.array([sys.argv[n] for n in range(3, 38, 2)])
reduced_seg_ids = [4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6]
reduced_locs = np.array([35.,35.,35.,35.,35.,35., 105.,105.,105.,105.,105.,105., 175.,175.,175.,175.,175.,175.])

data = {'reduced':[], 'detailed':[]}

fig, ax = plt.subplots()
trace_fig, trace_ax = plt.subplots()

for distance_index in range(n_points):
    seg_id_reduced = reduced_seg_ids[distance_index]
    seg_id_detailed = detailed_seg_ids[distance_index]
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
colors = ['blue', 'green']
time = np.loadtxt('{0}/time.dat'.format(sim_path))
maxima = {'reduced':{}, 'detailed':{}}
baselines = {'reduced':{}, 'detailed':{}}
amplitudes = {'reduced':{}, 'detailed':{}}
attenuation = {}


for k,cell_type in enumerate(data.keys()):
    for distance_index in range(n_points):
	trace_ax.plot(time,
		      data[cell_type][distance_index][1],
		      color=colors[k])
    maxima[cell_type]['soma'] = np.array([data[cell_type][dist_index][1].max() for dist_index in range(n_points)])
    maxima[cell_type]['dend'] = np.array([data[cell_type][dist_index][2].max() for dist_index in range(n_points)])
    baselines[cell_type]['soma'] = np.array([data[cell_type][dist_index][1][-1] for dist_index in range(n_points)])
    baselines[cell_type]['dend'] = np.array([data[cell_type][dist_index][2][-1] for dist_index in range(n_points)])
    amplitudes[cell_type]['soma'] = maxima[cell_type]['soma'] - baselines[cell_type]['soma']
    amplitudes[cell_type]['dend'] = maxima[cell_type]['dend'] - baselines[cell_type]['dend']

    attenuation[cell_type] = np.concatenate((np.array([1.]),
					     amplitudes[cell_type]['soma']/amplitudes[cell_type]['dend']))
    positions = np.array([0.] + [float(data[cell_type][dist_index][0]) for dist_index in range(n_points)])

    l = curve_fit(attenuation_function,
		  positions,
		  attenuation[cell_type],
		  [initial_l])[0]
    ax.plot(positions,
	    attenuation[cell_type],
	    marker='o', color=colors[k],
	    label=cell_type)
    print(cell_type, l)
    values = np.array([attenuation_function(x, l) for x in np.arange(0,200,0.1)])
    ax.plot(np.arange(0,200,0.1), values, color='black', linestyle=':')

ax.legend(loc='best')



plt.show()
