#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Single cell rate I/O plotting script. Usage: python rate_IO_plot.py iopf1415788530.0
import sys
import glob
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import seaborn as sns

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)

batch_ref = sys.argv[1]
stim_rate_range = range(1, 440, 40)
sim_duration = 6000. # ms
transient = 1000. # ms

fig, ax = plt.subplots()

cell_types = ['Solinas', 'Vervaeke', 'reduced']
rates = {}
threshold = 'min20'

for cell_type in cell_types:
    rates[cell_type] = []
    for amplitude in stim_rate_range:
	filename = '../simulations/{0}_{1}/Golgi_{2}_0.SPIKES_{3}.spike'.format(batch_ref,
										amplitude,
										cell_type,
										threshold)
        spikes = np.loadtxt(filename)
        rate = 1000 * np.sum(spikes > transient)/(sim_duration - transient) # (Hz)
        rates[cell_type].append(rate)


    ax.plot(stim_rate_range, rates[cell_type], marker='o', label=cell_type)

ax.legend(loc='best')
ax.set_xlabel('Stimulation rate (Hz)')
ax.set_ylabel('Firing rate (Hz)')

ax.locator_params(axis='x', tight=False, nbins=5)
ax.locator_params(axis='y', tight=False, nbins=5)
ax.set_xlim((-10, 410))
plt.tight_layout()
fig.savefig("fig/rate_IO_{}.pdf".format(batch_ref[2:4]))
plt.show()
