#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import random
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

timestamp = sys.argv[1] # 1422638890.67

rewiring_p_range = [1., 1e-1, 1e-2, 1e-3, 1e-4]
n_models = len(rewiring_p_range)
n_trials = 3

sim_duration = 2000
n_cells = 720
n_excluded_cells = 360
n_cells_in_raster = 45

fig, ax = plt.subplots(figsize=(4,3),ncols=1,nrows=2, sharex=True)

for i, (rewiring_p, color) in enumerate(zip(rewiring_p_range, sns.color_palette())):
    print i
    distances = []
    for trial in range(n_trials):
        cells_for_raster = random.sample(range(n_cells), n_cells_in_raster)
        cells_for_distance = random.sample(range(n_cells), n_cells - n_excluded_cells)

        sim_ref = utils.desynchronisation_small_world(timestamp,
                                                      rewiring_p,
                                                      trial)
        sim_dir = '../simulations/' + sim_ref
        time = np.loadtxt(sim_dir + '/time.dat')
        spike_trains = []

        raster_index = 0 
        for cell in range(n_cells):
            spikes = np.loadtxt('{}/Golgi_network_reduced_large_{}.SPIKES_min20.spike'.format(sim_dir, cell))
            spike_trains.append(pyspike.add_auxiliary_spikes(spikes, sim_duration))
            if trial==0 and cell in cells_for_raster:
                ax[0].scatter(spikes,
                              np.zeros_like(spikes)+(raster_index+i*n_cells_in_raster),
                              marker='|',
                              s=2,
                              c=color)
                raster_index += 1
        distances.append(pyspike.spike_profile_multi([spike_trains[c] for c in cells_for_distance]))
    
    # average synchrony index across trials
    average_distance = distances[0]
    for distance in distances[1:]:
        average_distance.add(distance)
    average_distance.mul_scalar(1./n_trials)
    x, y = average_distance.get_plottable_data()
    ax[1].plot(x, 1-y, lw=2, c=color, label='{:g}'.format(rewiring_p))

ax[1].legend(loc='lower right', title='Rewiring probability')
ax[0].locator_params(tight=True, nbins=4)
ax[1].locator_params(axis='y', tight=False, nbins=4)
ax[0].set_yticks([45 * k for k in range(n_models)])
ax[0].set_ylim((45*n_models-1, 0))
ax[1].set_xticks([0, 310, 1000, 2000])
ax[0].set_ylabel('Cell number')
ax[1].set_xlabel('Time (ms)')
ax[1].set_ylabel('Synchrony index')
plt.tight_layout()
fig.savefig('desynchronisation_small_world.pdf')

        
                    
        
