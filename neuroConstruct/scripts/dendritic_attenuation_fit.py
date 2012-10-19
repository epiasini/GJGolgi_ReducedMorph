#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates and compares attenuation for the detailed and reduced
morphology models of the Golgi cell.

Attenuation is defined in the following way: stimulate a dendrite with
an aEPSP far away from the soma (~200um), and measure the voltage
response at various points along the path connecting the stimulation
point to the soma. The attenuation is the ratio between the peak of
the somatic and the dendritic response.

This definition of attenuation is the one that's used in figure 3C in
Vervaeke2012.
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
rec_segs_detailed = [sys.argv[n] for n in range(2, 38, 2)]
rec_segs_reduced = [6,5,4,0]

distances = {'reduced': np.array([0., 35., 105., 175.]),
	     'detailed': np.array([sys.argv[n] for n in range(3, 38, 2)], dtype=np.float)}

sim_ref = timestamp
sim_path = '../simulations/' + sim_ref

traces = {'reduced':np.array([np.loadtxt('{0}/Golgi_reduced_0.{1}.dat'.format(sim_path,
									      seg_id))
			      for seg_id in rec_segs_reduced]),

	  'detailed':np.array([np.loadtxt('{0}/Golgi_Vervaeke_0.{1}.dat'.format(sim_path,
										seg_id))
			       for seg_id in rec_segs_detailed])}
attenuation = {}
space_constants = {}
att_fig, att_ax = plt.subplots()
for cell_type in ['reduced', 'detailed']:
    maxima = traces[cell_type].max(axis=1)
    baselines = traces[cell_type][:,-1]
    amplitudes = maxima - baselines
    attenuation[cell_type] = amplitudes/amplitudes[0]
    space_constants[cell_type] = curve_fit(attenuation_function,
					   distances[cell_type],
					   attenuation[cell_type],
					   [initial_l])[0][0]
    line = att_ax.plot(distances[cell_type], attenuation[cell_type],
		       marker='o', label=cell_type)[0]
    color = line.get_color()
    values = np.array([attenuation_function(x, space_constants[cell_type])
		       for x in np.arange(0,200,0.1)])
    att_ax.plot(np.arange(0,200,0.1), values, color=color, linestyle=':',
		label=r'$\lambda = {0:.1f}$'.format(space_constants[cell_type]))

att_ax.legend(loc='best')

plt.show()
