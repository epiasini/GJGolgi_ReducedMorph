#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Network I/O plotting script.
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt

timestamp = sys.argv[1]
stim_rate_range = range(1, 440, 40)

inn_fig, inn_ax = plt.subplots()
ninn_fig, ninn_ax = plt.subplots()

cell_types = ['reduced', 'Vervaeke']
rates = {}
sigmas = {}
threshold = 'min20'

for cell_type in cell_types:
    rates[cell_type] = {'innervated':[], 'non_innervated':[]}
    sigmas[cell_type] = {'innervated':[], 'non_innervated':[]}
    for amplitude in stim_rate_range:
	innervated_rates = []
	non_innervated_rates = []
	for cell_n in range(45):
	    filename='../simulations/{0}_{1}/Golgi_network_{2}_{3}.SPIKES_{4}.spike'.format(timestamp,
											    amplitude,
											    cell_type,
											    cell_n,
											    threshold)
	    try:
		cell_rate = len(np.loadtxt(filename))
	    except TypeError:
		cell_rate = 0
	    if cell_n < 15:
		innervated_rates.append(cell_rate)
	    else:
		non_innervated_rates.append(cell_rate)
	innervated_rates = np.array(innervated_rates)
	non_innervated_rates = np.array(non_innervated_rates)
	rates[cell_type]['innervated'].append(innervated_rates.mean())
	rates[cell_type]['non_innervated'].append(non_innervated_rates.mean())
	sigmas[cell_type]['innervated'].append(innervated_rates.std())
	sigmas[cell_type]['non_innervated'].append(non_innervated_rates.std())
    inn_ax.errorbar(stim_rate_range, rates[cell_type]['innervated'],
		    yerr=sigmas[cell_type]['innervated'], label=cell_type, fmt='o-')
    ninn_ax.errorbar(stim_rate_range, rates[cell_type]['non_innervated'],
		    yerr=sigmas[cell_type]['non_innervated'], label=cell_type, fmt='o-')


inn_ax.legend(loc='best')
inn_ax.set_xlabel('stimulation rate (Hz)')
inn_ax.set_ylabel('firing rate (Hz)')
inn_fig.suptitle('Direct innervated GoCs')
ninn_ax.legend(loc='best')
ninn_ax.set_xlabel('stimulation rate (Hz)')
ninn_ax.set_ylabel('firing rate (Hz)')
ninn_fig.suptitle('Non innervated GoCs')
plt.show()
