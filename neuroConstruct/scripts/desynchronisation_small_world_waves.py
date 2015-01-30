#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make a video showing desynchronisation in a small-world-like GoC network.
"""
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import utils


timestamp = '1422638890.67'
data_file_name = 'desynch_voltage.npy'

rewiring_p = 1e-3

sim_duration = 2000
sim_step = 0.05
video_step = 1.
sim_to_video_scale = video_step/sim_step


n_cells = 720

sim_ref = utils.desynchronisation_small_world(timestamp,
                                              rewiring_p,
                                              trial=0)
sim_dir = '../simulations/' + sim_ref
time_points = np.loadtxt(sim_dir + '/time.dat')
g = nx.read_graphml('/home/ucbtepi/thesis/data/GoC_net_structures/graph_' + sim_ref + '.graphml')

try:
    v = np.load(data_file_name)
except IOError:
    v = np.zeros((time_points.size, n_cells))
    for cell in range(n_cells):
        print("loading cell {}/{}".format(cell, n_cells))
        voltage_trace = np.loadtxt('{}/Golgi_network_reduced_large_{}.dat'.format(sim_dir, cell))
        v[:,cell] = voltage_trace
    np.save(data_file_name, v)

# normalise voltage
v = (v - v.min()) / (v.max() - v.min())

fig, ax = plt.subplots(figsize=(5,3))
pos = nx.spectral_layout(g)
node_collection = nx.draw_networkx_nodes(g, pos, ax=ax, node_size=5+50*(v[0]**2), zorder=2, node_color=v[0], cmap='RdYlBu_r', linewidths=0)
edge_collection = nx.draw_networkx_edges(g, pos, ax=ax, zorder=1)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

Blues = plt.get_cmap('RdYlBu_r')

def animate(i):
    print("animating frame {}/{}".format(i, time_points.size/sim_to_video_scale))
    node_collection.set_sizes(5+50*v[i*sim_to_video_scale])
    node_collection.set_color([Blues(x) for x in v[i*sim_to_video_scale]])
    return node_collection,
print('saving animation')
anim = animation.FuncAnimation(fig, animate,
                               frames=int(np.floor(time_points.size/sim_to_video_scale)),
                               interval=20,
                               blit=True)


#anim.save('basic_animation.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
fig.savefig('sw_graph.pdf')
