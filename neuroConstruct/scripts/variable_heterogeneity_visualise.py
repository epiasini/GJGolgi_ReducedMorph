#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt
import pyspike

import utils

timestamp = '1422357603.88'#'1422041472.12'#sys.argv[1]
deg_mean_range = [9.25]#, 12.40]
clustering_range = [0.52]#, 0.64]
n_trials = 1

sim_duration = 2000
n_cells = 45

for deg_mean in deg_mean_range:
    for clustering in clustering_range:
        for trial in range(n_trials):
            sim_ref = utils.variable_heterogeneity(timestamp,
                                                   deg_mean,
                                                   clustering,
                                                   trial)
            sim_dir = '../simulations/' + sim_ref
            time = np.loadtxt(sim_dir + '/time.dat')
            fig, ax = plt.subplots(ncols=1,nrows=2, sharex=True)
            spike_trains = []
            for cell in range(n_cells):
                spikes = np.loadtxt('{}/golgi_group_{}_0.SPIKE_min20.spike'.format(sim_dir,
                                                                                   cell))
                ax[0].scatter(spikes, np.zeros_like(spikes)+cell, marker='|', c='r')
                spike_trains.append(pyspike.add_auxiliary_spikes(spikes, sim_duration))
            x, y = pyspike.spike_profile_multi(spike_trains).get_plottable_data()
            ax[1].plot(x, 1-y, lw=2, c='k')
            fig.suptitle('{} {} {}'.format(deg_mean, clustering, trial))
            fig.savefig('variable_heterogeneity_test_{}_{}_{}.pdf'.format(deg_mean, clustering, trial))

        
                    
        
