#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Coupling strength plotting script for a the Golgi network connected
## though gap junctions as defined in Vervaeke2010 or Vervaeke2012.

import sys
import numpy as np
from matplotlib import pyplot as plt

import utils

timestamp = sys.argv[1]


cell_types = ['reduced']
gj_conn_types = ['2010', '2012']
n_trials = 1
averages = {}
sigmas = {}

for gj_conn_type in gj_conn_types:
    averages[gj_conn_type] = {}
    sigmas[gj_conn_type] = {}
    fig, ax = plt.subplots()
    for cell_type in cell_types:
	averages[gj_conn_type][cell_type] = []
	sigmas[gj_conn_type][cell_type] = []
	responses = []
	for trial in range(n_trials):
	    responses.append(np.zeros(shape=(45,45), dtype=np.float))
	    for stim_cell in range(45):
		sim_ref = utils.cs_sim_ref(timestamp,
					   gj_conn_type,
					   stim_cell,
					   trial)
		for rec_cell in range(45):
		    filename = '../simulations/{0}/Golgi_network_{1}_TTX_{2}.dat'.format(sim_ref,
											 cell_type,
											 rec_cell)
		    try:
			responses[-1][stim_cell, rec_cell] = np.loadtxt(filename)[-1]
		    except IOError:
			print('Data file not found: {0}'.format(sim_ref))
			responses[-1][stim_cell, rec_cell] = np.NaN
	print responses[0].shape
	print responses[0]
	#responses = np.array(responses).mean(axis=0)
        #averages[gj_conn_type][cell_type].append(responses.mean())
	#sigmas[gj_conn_type][cell_type].append(responses.std())
	#baseline = averages[gj_conn_type][cell_type][4]
	#averages[gj_conn_type][cell_type] = averages[gj_conn_type][cell_type] - baseline
	plot = ax.imshow(responses[0], interpolation='none', cmap='coolwarm')
	fig.colorbar(plot)


#ax.legend(loc='best')
#ax.set_xlabel('stimulation amplitude (pA)')
#ax.set_ylabel('voltage response (mV)')
#fig.suptitle('Steady-state current-voltage relations')
plt.show()
