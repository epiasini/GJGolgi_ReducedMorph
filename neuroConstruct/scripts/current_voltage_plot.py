#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Current-voltage relation plotting script for a single cell embedded
## in a golgi network though gap junctions as defined in Vervaeke2010
## or Vervaeke2012.

import sys
import numpy as np
from matplotlib import pyplot as plt

import utils

timestamp = sys.argv[1]
stim_range = range(-200, 220, 50)

fig, ax = plt.subplots()

cell_types = ['reduced', 'Vervaeke']
gj_conn_types = ['2010', '2012']
n_trials = 10
averages = {}
sigmas = {}

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
	ax.errorbar(stim_range, averages[gj_conn_type][cell_type],
		    yerr=sigmas[gj_conn_type][cell_type], label=cell_type+gj_conn_type, fmt='o-')

ax.legend(loc='best')
ax.set_xlabel('stimulation amplitude (pA)')
ax.set_ylabel('voltage response (mV)')
fig.suptitle('Steady-state current-voltage relations')
plt.show()
