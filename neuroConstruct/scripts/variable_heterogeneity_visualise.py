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

timestamp = sys.argv[1]

mean_scaling_range = [1.]
variance_scaling_range = [0.125, 1., 8.]#[0.125, 0.25, 0.5, 1., 2., 4., 8.]
n_trials = 1

sim_duration = 2000
n_cells = 45
n_stimulated_cells = 10

fig, ax = plt.subplots(ncols=1,nrows=2, sharex=True)

for i, ((mean_scaling, variance_scaling), color) in enumerate(zip(itertools.product(mean_scaling_range, variance_scaling_range), sns.color_palette())):
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
            ax[0].scatter(spikes,
                          np.zeros_like(spikes)+cell+i*n_cells,
                          marker='|',
                          s=1,
                          c=color)
            spike_trains.append(pyspike.add_auxiliary_spikes(spikes, sim_duration))
        x, y = pyspike.spike_profile_multi(spike_trains[n_stimulated_cells:]).get_plottable_data()
        ax[1].plot(x, 1-y, lw=2, c=color, label='m {:g}, v {:g}'.format(mean_scaling,
                                                                        variance_scaling))

ax[1].legend(loc='best')
fig.savefig('variable_heterogeneity.pdf')

        
                    
        
