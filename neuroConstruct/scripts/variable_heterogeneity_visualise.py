#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt

import utils

timestamp = sys.argv[1]
deg_mean_range = [35]
deg_sigma_cv_range = [1./3]
syn_strength_noise = 0.6

for deg_mean in deg_mean_range:
    deg_sigma_range = [x*deg_mean for x in deg_sigma_cv_range]
    for deg_sigma in deg_sigma_range:
	sim_ref = utils.variable_heterogeneity(timestamp,
					       deg_mean,
					       deg_sigma)
	sim_dir = '../simulations/' + sim_ref
	time = np.loadtxt(sim_dir + '/time.dat')
	fig, ax = plt.subplots()
	for datafile in glob.glob(sim_dir + '/Golgi_*.dat'):
	    ax.plot(time, np.loadtxt(datafile))
	fig.suptitle('{0} {1} {2}'.format(deg_mean, deg_sigma, syn_strength_noise))

plt.show()
