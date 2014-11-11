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
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import seaborn as sns

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{euler}']
matplotlib.rcParams['font.sans-serif'].insert(0, 'Bitstream Vera Sans')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.size'] = 8.0
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['xtick.labelsize'] = 'medium'
matplotlib.rcParams['ytick.labelsize'] = 'medium'
matplotlib.rcParams['axes.labelsize'] = 'medium'
figsize=(3.5,3)

def attenuation_function(x, l):
    return np.exp(-x/l)

initial_a = 1.
initial_l = 46.

n_points = 18
sim_ref = sys.argv[1]
rec_segs_detailed = [sys.argv[n] for n in range(2, 38, 2)]
rec_segs_reduced = [6,5,4,0]

distances = {'reduced': np.array([0., 35., 105., 175.]),
	     'Vervaeke': np.array([sys.argv[n] for n in range(3, 38, 2)], dtype=np.float)}

sim_path = '../simulations/' + sim_ref

traces = {'reduced':np.array([np.loadtxt('{0}/Golgi_reduced_0.{1}.dat'.format(sim_path,
									      seg_id))
			      for seg_id in rec_segs_reduced]),

	  'Vervaeke':np.array([np.loadtxt('{0}/Golgi_Vervaeke_0.{1}.dat'.format(sim_path,
										seg_id))
			       for seg_id in rec_segs_detailed])}
attenuation = {}
space_constants = {}
fig, ax = plt.subplots(figsize=figsize)
for cell_type in ['reduced', 'Vervaeke']:
    maxima = traces[cell_type].max(axis=1)
    baselines = traces[cell_type][:,-1]
    amplitudes = maxima - baselines
    attenuation[cell_type] = amplitudes/amplitudes[0]
    space_constants[cell_type] = curve_fit(attenuation_function,
					   distances[cell_type],
					   attenuation[cell_type],
					   [initial_l])[0][0]
    line = ax.plot(distances[cell_type], attenuation[cell_type],
		       marker='o', label=cell_type, linestyle='')[0]
    color = line.get_color()
    values = np.array([attenuation_function(x, space_constants[cell_type])
		       for x in np.arange(0,200,0.1)])
    ax.plot(np.arange(0,200,0.1), values, color=color, linestyle='-',
		label=r'$\lambda = {0:.1f}$'.format(space_constants[cell_type]))

ax.set_xlabel(r"Distance ($\mu$m)")
ax.set_ylabel("Attenuation (a.u.)")
ax.legend(loc='best')
ax.locator_params(tight=True, nbins=5)
plt.tight_layout()

fig.savefig("fig/dendritic_attenuation.pdf")

plt.show()
