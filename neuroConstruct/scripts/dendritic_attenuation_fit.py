#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates and compares attenuation for the detailed and reduced
morphology models of the Golgi cell.

Attenuation is defined in the following way: stimulate a dendrite with
an aEPSP far away from the soma (~200um), and measure the voltage
response at various points along the path connecting the stimulation
point to the soma. The attenuation is the ratio between the peak of
the dendritic and the somatic response.

This definition of attenuation is the one that's used in figure 3C in
Vervaeke2012.

Usage: pythonv dendritic_attenuation_fit.py da1415731151.61 564 0.0 549 10.5919850469 540 19.8728369474 506 31.4311226606 493 42.2954727411 481 48.7114700973 447 60.1053161472 435 70.3048134446 424 78.4546051621 406 92.3413149714 378 100.099774525 360 111.302530318 348 118.78751111 328 131.306678474 315 138.444323123 205 149.846530199 188 159.575384408 139 171.218757302 124 180.70650842 0 208.072402954
"""
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import seaborn as sns

rc = matplotlib.rc_params_from_file('/home/ucbtepi/thesis/matplotlibrc.thesis',
                                    use_default_template=False)
matplotlib.rcParams.update(rc)

def attenuation_function(x, l):
    return np.exp(-x/l)

initial_a = 1.
initial_l = 46.

n_points = 18
sim_ref = sys.argv[1]
rec_segs_detailed = [sys.argv[n] for n in range(2, 38, 2)]
rec_segs_reduced = [6,5,4,0]

distances = {'reduced': np.array([0., 35., 105., 175.]),
	     'Vervaeke': np.array([sys.argv[n] for n in range(3, 38, 2)], dtype=np.float)}

sim_path = '../simulations/' + sim_ref

traces = {'reduced':np.array([np.loadtxt('{0}/Golgi_reduced_0.{1}.dat'.format(sim_path,
									      seg_id))
			      for seg_id in rec_segs_reduced]),

	  'Vervaeke':np.array([np.loadtxt('{0}/Golgi_Vervaeke_0.{1}.dat'.format(sim_path,
										seg_id))
			       for seg_id in rec_segs_detailed])}
attenuation = {}
space_constants = {}
fig, ax = plt.subplots()
for cell_type in ['reduced', 'Vervaeke']:
    maxima = traces[cell_type].max(axis=1)
    baselines = traces[cell_type][:,-1]
    amplitudes = maxima - baselines
    attenuation[cell_type] = amplitudes/amplitudes[0]
    space_constants[cell_type] = curve_fit(attenuation_function,
					   distances[cell_type],
					   attenuation[cell_type],
					   [initial_l])[0][0]
    line = ax.plot(distances[cell_type], attenuation[cell_type],
		       marker='o', label=cell_type, linestyle='')[0]
    color = line.get_color()
    values = np.array([attenuation_function(x, space_constants[cell_type])
		       for x in np.arange(0,200,0.1)])
    ax.plot(np.arange(0,200,0.1), values, color=color, linestyle='-',
		label=r'$\lambda = {0:.1f}$'.format(space_constants[cell_type]))

ax.set_xlabel(r"Distance ($\mu$m)")
ax.set_ylabel("Attenuation (a.u.)")
ax.legend(loc='best')
ax.locator_params(tight=False, nbins=5)
ax.set_xlim((-5, 200))
ax.set_ylim((0, 1.05))
plt.tight_layout()

fig.savefig("fig/dendritic_attenuation.pdf")

plt.show()
