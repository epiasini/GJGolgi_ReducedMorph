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
rewiring_p_range = [1., 0.1, 0.01, 0.001]
make_animation = True


sim_duration = 2000
sim_step = 0.05
video_step = 0.5
sim_to_video_scale = video_step/sim_step


n_cells = 720


for rewiring_p in rewiring_p_range:
    sim_ref = utils.desynchronisation_small_world(timestamp,
                                                  rewiring_p,
                                                  trial=0)
    sim_dir = '/Users/eugenio/mnt/virtual/thesis_edit'
    time_points = np.loadtxt(sim_dir + '/time.dat')
    graphml_file_name = sim_dir + '/graph_' + sim_ref + '.graphml'
    data_file_name = sim_dir + '/desynch_voltage_{}.npy'.format(sim_ref)


    g = nx.read_graphml(graphml_file_name)

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
    color = (v - v.min()) / (-30. - v.min())
    color[color>1] = 1
    size = (v - v.min()) / (-10. - v.min())
    size[size>1] = 1
    size *= 50

    if rewiring_p==1:
        pos = nx.random_layout(g)
    else:
        pos = nx.spectral_layout(g)
    if make_animation:
        fig, ax = plt.subplots(figsize=(5,3))
        node_collection = nx.draw_networkx_nodes(g, pos, ax=ax, node_size=size[0], zorder=2, node_color=color[0], cmap='RdYlBu_r', linewidths=0)
        edge_collection = nx.draw_networkx_edges(g, pos, ax=ax, zorder=1, width=0.75, alpha=0.5, edge_color='k')
    else:
        fig, ax = plt.subplots(figsize=(2,2.5))
        node_collection = nx.draw_networkx_nodes(g, pos, ax=ax, node_size=10, zorder=2, node_color='#7fbc41', linewidths=0.3, alpha=1)
        node_collection.set_edgecolor('#4d9221')
        edge_collection = nx.draw_networkx_edges(g, pos, ax=ax, zorder=1, width=0.75, alpha=0.5, edge_color='k')#'#bababa')

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.locator_params(tight=True, nbins=4)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

    pts = np.array([[0.88,0], [1,0], [1,0.2]])
    polygon = plt.Polygon(pts, closed=True, fill=True, linewidth=0, color='#feb24c')
    #circle = plt.Circle((0.95,0.05),.05,color='#feb24c')
    #fig.gca().add_artist(circle1)


    Blues = plt.get_cmap('RdYlBu_r')

    if make_animation:
        def animate(i):
            print("animating frame {}/{}".format(i, time_points.size/sim_to_video_scale))
            node_collection.set_sizes(size[i*sim_to_video_scale])
            node_collection.set_color([Blues(c) for c in color[i*sim_to_video_scale]])
            result = node_collection,
            if rewiring_p==0.001:
                if 310 <= i*video_step <= 312 or 315 <= i*video_step <= 317:

                    #circle = plt.Circle((0.8,0.2),.2,color='#feb24c')
                    if polygon in ax.patches:
                        print ax.patches
                        polygon.remove()
                    ax.add_patch(polygon)
                    result = node_collection, polygon
                elif (312 < i*video_step < 315 or 317 < i*video_step):
                    if polygon in ax.patches:
                        polygon.remove()
                    result = node_collection, polygon
            return result
        print('saving animation')
        anim = animation.FuncAnimation(fig, animate,
                                       frames=int(np.floor(time_points.size/sim_to_video_scale)),
                                       interval=20,
                                       blit=True)


        anim.save('desynchronisation_{}.mp4'.format(sim_ref), fps=40, extra_args=['-vcodec', 'libx264'], dpi=200)
    fig.savefig('sw_graph_{}.pdf'.format(sim_ref))
