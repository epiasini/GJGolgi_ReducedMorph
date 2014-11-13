#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Network I/O plotting script, looking at basolateral or apical
## stimulation at a subset of cells, and how this propagates through
## gap junctions to cells that are not directly innervated.
## usage: pythonv rate_IO_network_plot.py net_1415883300.46
import sys
import numpy as np

import matplotlib
import seaborn as sns
rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)

figsize = (2.5, 3)

batch_ref = sys.argv[1]
stim_rate_range = range(1, 440, 40)
n_points = len(stim_rate_range)

fig, ax = matplotlib.pyplot.subplots((1,2), figsize=figsize)

cell_types = ['Vervaeke', 'reduced']
colors = {'Vervaeke': '#55A868', 'reduced': '#4C72B0'}
threshold = 'min20'

zorder_counter = 1
lines = []

for stim_source, lin zip(['pf'], ['-']):
    ax.set_color_cycle(None) # reset color cycle
    for cell_type in cell_types:
        color = colors[cell_type]
        rates = {'innervated': np.zeros(n_points), 'non_innervated': np.zeros(n_points)}
        sigmas = {'innervated': np.zeros(n_points), 'non_innervated': np.zeros(n_points)}
        for i, amplitude in enumerate(stim_rate_range):
            sim_ref = 'net_{}{}_{}'.format(stim_source,
                                           batch_ref.lstrip('net_'),
                                           amplitude)
            innervated_rates = []
            non_innervated_rates = []
            for cell_n in range(45):
                filename='../simulations/{0}/Golgi_network_{1}_{2}.SPIKES_{3}.spike'.format(sim_ref, cell_type, cell_n, threshold)
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
            rates['innervated'][i] = innervated_rates.mean()
            rates['non_innervated'][i] = non_innervated_rates.mean()
            sigmas['innervated'][i] = innervated_rates.std()
            sigmas['non_innervated'][i] = non_innervated_rates.std()
        ax.plot(stim_rate_range, rates['innervated'], label="{}".format(cell_type), linestyle=linestyle, linewidth=1.5, color=color)
        ax.fill_between(stim_rate_range, rates['innervated']-sigmas['innervated'], rates['innervated']+sigmas['innervated'], color=color, alpha=0.3)
        ax.plot(stim_rate_range, rates['non_innervated'], linestyle="--", linewidth=1.5, dashes=(4,1.5), color=color)
        ax.fill_between(stim_rate_range, rates['non_innervated']-sigmas['non_innervated'], rates['non_innervated']+sigmas['non_innervated'], color=color, alpha=0.3)
        # ax.errorbar(stim_rate_range, rates[cell_type]['innervated'],
        #             yerr=sigmas[cell_type]['innervated'],
        #             label="{}, innervated".format(cell_type), fmt=linestyle, zorder=zorder_counter, capsize=1.5, capthick=0.5, elinewidth=1, alpha=0.7)
        # ax.errorbar(stim_rate_range, rates[cell_type]['non_innervated'],
        #             yerr=sigmas[cell_type]['non_innervated'],
        #             label="{}, non innervated".format(cell_type), fmt=linestyle, zorder=zorder_counter+1, capsize=1.5, capthick=0.5, elinewidth=1, alpha=0.7)
        zorder_counter += 1


ax.legend(loc='center right')
ax.set_xlabel('Stimulation rate (Hz)')
ax.set_ylabel('Firing rate (Hz)')
#fig.suptitle('Direct innervated GoCs')
ax.locator_params(axis='x', tight=False, nbins=4)
ax.locator_params(axis='y', tight=False, nbins=4)
ax.set_xlim((-10, 410))
ax.set_ylim((0, 80))
matplotlib.pyplot.tight_layout()

if not 10 in ax.get_yticks():
    ax.set_yticks(np.hstack((ax.get_yticks(), 10)))

fig_filename = 'fig/rate_IO_network.pdf'
fig.savefig(fig_filename)
print("result saved to {}".format(fig_filename))
matplotlib.pyplot.show()
