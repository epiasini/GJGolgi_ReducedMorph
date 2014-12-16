#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Current-voltage relation plotting script for a single cell embedded
## in a golgi network though gap junctions as defined in Vervaeke2010
## or Vervaeke2012. Usage: pythonv current_voltage_plot.py 1415727840.3

import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import seaborn as sns

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)
figsize = (3.5, 2)

import utils

timestamp = sys.argv[1]
stim_range = range(-200, 220, 50)

fig, ax = plt.subplots(figsize=figsize)

cell_types = ['Vervaeke','reduced']
gj_conn_types = ['2010', '2012']
colors = ['#55A868', '#4C72B0']
dash_styles = [(4,1.5), (None, None)]
n_trials = 1
averages = {}
sigmas = {}
lines = []

for gj_conn_type, dashes in zip(gj_conn_types, dash_styles):
    averages[gj_conn_type] = {}
    sigmas[gj_conn_type] = {}
    for cell_type, color in zip(cell_types, colors):
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
            lines.extend(ax.plot(stim_range, averages[gj_conn_type][cell_type], marker='o', c=color, dashes=dashes))
        
lines[2].set_label("Vervaeke")
lines[3].set_label("reduced")
ax.legend(loc='lower right')
ax.set_xlabel('Injected current (pA)')
ax.set_ylabel('Voltage response (mV)')
#fig.suptitle('Steady-state current-voltage relations')
ax.locator_params(tight=False, nbins=5)
ax.set_xlim((-210, 210))
plt.tight_layout()
fig.savefig("fig/current_voltage.pdf")
plt.show()
