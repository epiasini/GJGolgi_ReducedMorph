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
matplotlib.rcParams.update(rc)

def fermi_function(r, a, r_0, delta):
    return a/(np.exp((r-r_0)/delta) + 1)

def vervaeke2010_conn_prob(r):
    result = -17.45 + 18.36 / (np.exp((r - 267.)/39.) + 1)
    result[result<0] = 0
    return result


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
print("Parameters for best fit Fermi function connection probability model:\n  a: {}\n  r₀: {}\n  Δ: {}".format(a, r_0, delta))

# plot
fig, ax = plt.subplots(figsize=(2.5,2))
x_values = np.linspace(0, 200, 1000)
# vervaeke 2010 model
ax.plot(x_values, vervaeke2010_conn_prob(x_values), color=seaborn.color_palette()[1], alpha=1, linewidth=2)
# fermi function model
ax.plot(x_values, fermi_function(x_values, a, r_0, delta), color=seaborn.color_palette()[2], alpha=1, linewidth=2)
# experiment
ax.plot(bin_centers.squeeze(), connection_probabilities.squeeze(), alpha=0.8, linewidth=1, ls='', marker='o')

ax.set_xlabel(r"Distance (\si{\micro\metre})")
ax.set_ylabel("Connection probability")
ax.locator_params(tight=False, nbins=3)
ax.set_xticks(np.insert(ax.get_xticks(), 1, 30))
plt.tight_layout()

fig.savefig("fig/GoC_GJ_connection_probability_fermi_fit.pdf")

