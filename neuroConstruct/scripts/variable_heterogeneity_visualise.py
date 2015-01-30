#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
import pyspike
import itertools

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)

import utils

timestamp = sys.argv[1] # 1422612003.0

mean_scaling_range = [1.]
variance_scaling_range = [1e-6, 1, 4, 9]
n_models = len(variance_scaling_range)
n_trials = 32

sim_duration = 2000
n_cells = 45
n_excluded_cells = 2

fig, ax = plt.subplots(figsize=(4,3),ncols=1,nrows=2, sharex=True)


for i, ((mean_scaling, variance_scaling), color) in enumerate(zip(itertools.product(mean_scaling_range, variance_scaling_range), sns.color_palette())):
    distances = []
    for trial in range(n_trials):
        sim_ref = utils.variable_heterogeneity(timestamp,
                                               mean_scaling,
                                               variance_scaling,
                                               trial)
        sim_dir = '../simulations/' + sim_ref
        time = np.loadtxt(sim_dir + '/time.dat')
        spike_trains = []
        for cell in range(n_cells):
            spikes = np.loadtxt('{}/Golgi_network_reduced_{}.SPIKES_min20.spike'.format(sim_dir, cell))
            spike_trains.append(pyspike.add_auxiliary_spikes(spikes, sim_duration))
            if trial==3:
                ax[0].scatter(spikes,
                              np.zeros_like(spikes)+(cell+i*n_cells),
                              marker='|',
                              s=2,
                              c=color)
        distances.append(pyspike.spike_profile_multi(spike_trains[n_excluded_cells:]))

    # average synchrony index across trials
    average_distance = distances[0]
    for distance in distances[1:]:
        average_distance.add(distance)
    average_distance.mul_scalar(1./n_trials)
    x, y = average_distance.get_plottable_data()
    ax[1].plot(x, 1-y, lw=2, c=color, label='{:g}'.format(variance_scaling))

ax[1].legend(loc='lower right', title='Variance scaling')
ax[0].locator_params(tight=True, nbins=4)
ax[1].locator_params(axis='y', tight=False, nbins=4)
ax[0].set_yticks([45 * k for k in range(n_models)])
ax[0].set_ylim((45*n_models-1, 0))
ax[1].set_xticks([0, 310, 1000, 2000])
ax[0].set_ylabel('Cell number')
ax[1].set_xlabel('Time (ms)')
ax[1].set_ylabel('Synchrony index')
plt.tight_layout()
fig.savefig('variable_heterogeneity.pdf')

        
                    
        
