#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Coupling strength plotting script for a the Golgi network connected
## though gap junctions as defined in Vervaeke2010 or Vervaeke2012.

import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import seaborn

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)

def fermi_function(x, a, r_0, delta):
    return a/(np.exp((x-r_0)/delta) + 1)


# load experimental data
exp_distances = np.loadtxt('../dataSets/koen_data/golgi_pair_distances.csv')
exp_couplings = np.loadtxt('../dataSets/koen_data/golgi_pair_couplings.csv')/100.
average_distances = np.mean(np.vstack([np.atleast_2d(exp_distances[0::2]),
                                       np.atleast_2d(exp_distances[1::2])]),
                            axis=0)
average_couplings = np.mean(np.vstack([np.atleast_2d(exp_couplings[0::2]),
                                       np.atleast_2d(exp_couplings[1::2])]),
                            axis=0)

# compute estimate of connection probability as a function of cell distance
total_pairs = np.histogram(average_distances, range=[exp_distances.min(),
                                                     exp_distances.max()])
connected_pairs = np.histogram(average_distances[average_couplings>0.01],
                               range=[exp_distances.min(),
                                      exp_distances.max()])
# ensure the histograms have the same bin edges
assert (total_pairs[1] == connected_pairs[1]).all()
bin_edges = total_pairs[1].reshape(total_pairs[1].size, -1)
bin_centers = np.mean(np.hstack([bin_edges[0:-1], bin_edges[1:]]),
                      axis=1).reshape(total_pairs[0].size, -1)
connection_probabilities = (connected_pairs[0].astype(np.float)/total_pairs[0].astype(np.float)).reshape(total_pairs[0].size,-1)

# fit Fermi function model
a, r_0, delta = curve_fit(fermi_function,
                          bin_centers.transpose().squeeze(),
                          connection_probabilities.transpose().squeeze(),
                          p0=[0.856, 122., 16.9])[0]

# plot
fig, ax = plt.subplots(figsize=(3,2.5))
ax.fill_between(bin_centers.squeeze(), connection_probabilities.squeeze(), alpha=0.7, linewidth=2)
x_values = np.linspace(bin_centers[0], 200, 1000)
ax.fill_between(x_values, fermi_function(x_values, a, r_0, delta), color=seaborn.color_palette()[2], alpha=0.6, linewidth=2)

ax.set_xlabel("Distance ($\mu$m)")

ax.locator_params(tight=False, nbins=3)
plt.tight_layout()

fig.savefig("fig/connection_probability.pdf")

