#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt
import pyspike

import utils

timestamp = sys.argv[1]

sim_duration = 2000
n_cells = 45
n_stimulated_cells = 10
years = '2010', '2012'
colors = 'r', 'b'

fig, ax = plt.subplots(ncols=1,nrows=2, sharex=True)


for k, (year, color) in enumerate(zip(years, colors)):

    sim_ref = 'rd_' + timestamp + '_y' + year
    sim_dir = '../simulations/' + sim_ref
    time = np.loadtxt(sim_dir + '/time.dat')
    spike_trains = []
    for cell in range(n_cells):
        spikes = np.loadtxt('{}/Golgi_network_reduced_{}.SPIKES_min20.spike'.format(sim_dir,
                                                                                    cell))
        ax[0].scatter(spikes, np.zeros_like(spikes)+cell+k*n_cells, marker='|', c=color)
        spike_trains.append(pyspike.add_auxiliary_spikes(spikes, sim_duration))
    x, y = pyspike.spike_profile_multi(spike_trains[n_stimulated_cells:]).get_plottable_data()
    ax[1].plot(x, 1-y, lw=2, c=color, label=year)

ax[1].legend(loc='best')
fig.savefig('reduced_desynchronisation.pdf')

        
                    
        
