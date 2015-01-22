#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt

import utils

timestamp = '1421889230.82'#sys.argv[1]
deg_mean_range = [15]
n_trials = 1

n_cells = 45

for deg_mean in deg_mean_range:
    for trial in range(n_trials):
	sim_ref = utils.variable_heterogeneity(timestamp,
					       deg_mean,
					       trial)
	sim_dir = '../simulations/' + sim_ref
	time = np.loadtxt(sim_dir + '/time.dat')
	fig, ax = plt.subplots()
        for cell in range(n_cells):
            spikes = np.loadtxt('{}/golgi_group_{}_0.SPIKE_min20.spike'.format(sim_dir,
                                                                               cell))
            ax.scatter(spikes, np.zeros_like(spikes)+cell, marker='|', c='r')
	fig.suptitle('{0} {1}'.format(deg_mean, trial))
        fig.savefig('variable_heterogeneity_test_{}_{}.pdf'.format(deg_mean, trial))
