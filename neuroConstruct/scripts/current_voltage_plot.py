#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Current-voltage relation plotting script for a single cell embedded
## in a golgi network though gap junctions as defined in Vervaeke2010
## or Vervaeke2012.

import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

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

import utils

timestamp = sys.argv[1]
stim_range = range(-200, 220, 50)

fig, ax = plt.subplots(figsize=figsize)
colors = sns.xkcd_palette(["windows blue", "amber", "faded green", "dusty purple"])
#sns.set_palette("colorblind")
#sns.palplot(sns.xkcd_palette(colors))
#colors = sns.palplot()

cell_types = ['Vervaeke','reduced']
gj_conn_types = ['2010', '2012']
n_trials = 1
averages = {}
sigmas = {}

color_index = 0

for gj_conn_type in gj_conn_types:
    averages[gj_conn_type] = {}
    sigmas[gj_conn_type] = {}
    for cell_type in cell_types:
	averages[gj_conn_type][cell_type] = []
	sigmas[gj_conn_type][cell_type] = []
	for stim_level in stim_range:
	    responses = []
	    for trial in range(n_trials):
		sim_ref = utils.ir_sim_ref(timestamp,
					   gj_conn_type,
					   stim_level,
					   trial)
		filename = '../simulations/{0}/Golgi_network_{1}_TTX_0.dat'.format(sim_ref,
										   cell_type)
		try:
		    responses.append(np.loadtxt(filename)[-1])
		except IOError:
		    print('Data file not found: {0}'.format(sim_ref))
	    responses = np.array(responses)
	    averages[gj_conn_type][cell_type].append(responses.mean())
	    sigmas[gj_conn_type][cell_type].append(responses.std())
	baseline = averages[gj_conn_type][cell_type][4]
	averages[gj_conn_type][cell_type] = averages[gj_conn_type][cell_type] - baseline
        label = '{}, {} net'.format(cell_type, gj_conn_type)
        if n_trials > 1:
            ax.errorbar(stim_range, averages[gj_conn_type][cell_type],
                        yerr=sigmas[gj_conn_type][cell_type], label=label, fmt='o-')
        else:
            ax.plot(stim_range, averages[gj_conn_type][cell_type], label=label,marker='o',linestyle='-', c=colors[color_index])
            color_index += 1
        

ax.legend(loc='best')
ax.set_xlabel('Stimulation amplitude (pA)')
ax.set_ylabel('Voltage response (mV)')
#fig.suptitle('Steady-state current-voltage relations')
ax.locator_params(tight=True, nbins=5)
plt.tight_layout()
fig.savefig("fig/current_voltage.pdf")
plt.show()
