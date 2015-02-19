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
import networkx as nx

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)

import utils

timestamp = '1422638890.67' # sys.argv[1]

rewiring_p_range = [1., 1e-1, 1e-2, 1e-3, 1e-4]
n_models = len(rewiring_p_range)
n_trials = 3

sim_duration = 2000
n_cells = 720
n_excluded_cells = 0
n_cells_in_raster = 45

fig, ax = plt.subplots(figsize=(4,3),ncols=1,nrows=2, sharex=True)
lines = []
labels = []

synchrony_indices = []

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
    synchrony_indices.append(1-average_distance.avrg())
    xmin = 50
    xmax = 1900
    x, y = average_distance.get_plottable_data()
    ximin = np.searchsorted(x, xmin)
    ximax = np.searchsorted(x, xmax)
    lines.append(ax[1].plot(x[ximin:ximax+1], 1-y[ximin:ximax+1], lw=2, c=color)[0])
    labels.append('{:g}'.format(rewiring_p))
    ax[1].plot(x[:ximin+1], 1-y[:ximin+1], lw=2, c=color, alpha=0.4)
    ax[1].plot(x[ximax:], 1-y[ximax:], lw=2, c=color, alpha=0.4)


#ax[1].legend(loc='lower right', title='Rewiring probability')
ax[0].locator_params(tight=True, nbins=4)
ax[1].locator_params(axis='y', tight=False, nbins=4)
ax[0].set_yticks([45 * k for k in range(n_models)])
ax[0].set_ylim((45*n_models-1, 0))
ax[1].set_xticks([0, 310, 1000, 2000])
ax[0].set_ylabel('Cell number')
ax[1].set_xlabel('Time (ms)')
ax[1].set_ylabel('Synchrony index')
plt.tight_layout()
fig.subplots_adjust(hspace=.5)
fig.legend(lines, labels, title='Rewiring probability', loc='center', ncol=n_models,
           bbox_to_anchor=(0.55, 0.55))
fig.savefig('desynchronisation_small_world.pdf')


# plot average clustering, mean path length and synchrony index as
# functions of the rewiring probability.

fig, ax = plt.subplots(figsize=(3, 2.5))
ax.set_xscale('log')
ax2 = ax.twinx()

clustering = np.zeros_like(rewiring_p_range)
path_length = np.zeros_like(rewiring_p_range)
for k, p in enumerate(rewiring_p_range):
    print p
    c = []
    l = []
    for t in range(5):
        g = nx.connected_watts_strogatz_graph(n_cells, 15, p, tries=100)
        c.append(np.array(nx.clustering(g).values()).mean())
        l.append(nx.average_shortest_path_length(g))
    clustering[k] = np.array(c).mean()
    path_length[k] = np.array(l).mean()

# compute reference clustering and path length values
c = []
l = []
for t in range(5):
    g0 = nx.connected_watts_strogatz_graph(n_cells, 15, 0, tries=100)
    c.append(np.array(nx.clustering(g0).values()).mean())
    l.append(nx.average_shortest_path_length(g0))
c0 = np.array(c).mean()
l0 = np.array(l).mean()    

ax.plot(rewiring_p_range, clustering/c0, c=sns.color_palette()[0], lw=2, marker='o', ls='--', zorder=10)
ax.plot(rewiring_p_range, path_length/l0, c=sns.color_palette()[0], lw=2, marker='o', zorder=11)
ax2.plot(rewiring_p_range, synchrony_indices, c=sns.color_palette()[2], lw=2, marker='o', zorder=12)

ax.set_ylim((0, 1.05))
ax2.set_ylim(0.8, 1.01)

ax.set_yticks([0., 0.5, 1.])
ax2.set_yticks([0.8, 0.9, 1.])

for tl in ax.get_yticklabels():
    tl.set_color(sns.color_palette()[0])
for tl in ax2.get_yticklabels():
    tl.set_color(sns.color_palette()[2])

ax.set_xlabel('Rewiring probability')
plt.tight_layout()
fig.savefig('desynchronisation_small_world_summary.pdf')
