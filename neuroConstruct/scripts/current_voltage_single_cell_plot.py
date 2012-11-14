#! /usr/bin/env python
# -*- coding: utf-8 -*- Current-voltage relation plotting script for a
## single Golgi cell removed from the network (in particolar, without
## gap junctions).

import sys
import numpy as np
from matplotlib import pyplot as plt

import utils

timestamp = sys.argv[1]
stim_range = range(-200, 220, 50)

fig, ax = plt.subplots()

cell_types = ['reduced', 'Vervaeke']
voltages = {}

for cell_type in cell_types:
    voltages[cell_type] = np.zeros(shape=len(stim_range), dtype=np.float)
    for k,stim_level in enumerate(stim_range):
	sim_ref = utils.ir_single_cell_sim_ref(timestamp,
					       stim_level)
	filename = '../simulations/{0}/Golgi_{1}_TTX_0.dat'.format(sim_ref,
								   cell_type)
	try:
	    voltages[cell_type][k] = np.loadtxt(filename)[-1]
	except IOError:
	    print('Data file not found: {0}'.format(sim_ref))

    baseline = voltages[cell_type][4]
    voltages[cell_type] = voltages[cell_type] - baseline
    ax.plot(stim_range, voltages[cell_type], label=cell_type, marker='o')

ax.legend(loc='best')
ax.set_xlabel('stimulation amplitude (pA)')
ax.set_ylabel('voltage response (mV)')
fig.suptitle('Steady-state current-voltage relations, single cell (no gap junctions).\ninput resistance (-200pA to 0pA): {0}M$\Omega$'.format(1000*voltages['reduced'][0]/float(stim_range[0])))
plt.show()
