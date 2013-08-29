#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Coupling strength plotting script for a the Golgi network connected
## though gap junctions as defined in Vervaeke2010 or Vervaeke2012.

import sys
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from scipy.spatial import distance

import utils

def entropy(p):
    return -(p[p>0]*np.log2(p[p>0])).sum()

def load_ordered_edge_list(timestamp, gj_conn_type, trial):
    edge_list = np.loadtxt(utils.cs_edge_list_file(timestamp, gj_conn_type, trial),
                           delimiter=",",
                           dtype=np.int)
    edge_list.sort()
    return edge_list

def load_unique_edge_list(timestamp, gj_conn_type, trial):
    edge_list = load_ordered_edge_list(timestamp, gj_conn_type, trial)
    return np.array(list(set([tuple(x) for x in edge_list])))

def load_edge_lengths(timestamp, gj_conn_type, trial):
    positions = np.loadtxt(utils.cs_cell_positions_file(timestamp, gj_conn_type, trial),
                           delimiter=",")
    distances = distance.squareform(distance.pdist(positions))
    return distances

# basic controls
timestamp = sys.argv[1]
gj_conn_types = ['2010', '2012']
n_cells = 45
n_trials = 1

# initialise data structures for storing results
coupling_coefficients = {} # coupling coefficients of all cell pairs
cc_conn = {} # coupling coefficients of anatomically connected cell pairs
edge_lengths = {} # pairwise distances between all cell pairs
graphs = {} # abstract graphs representing network realisations
degrees = {} # degree sequences

# plotting objects
colours = {'2010':'k', '2012':'g'}
cc_axes = {}
cc_ims = {}
cc_fig, (cc_axes['2010'], cc_axes['2012']) = plt.subplots(ncols=2)
cc_vs_d_fig, cc_vs_d_ax = plt.subplots()
deg_hist_fig, deg_hist_ax = plt.subplots()
cc_hist_fig, cc_hist_ax = plt.subplots()

for gj_conn_type in gj_conn_types:
    responses = []
    graphs[gj_conn_type] = []
    degrees[gj_conn_type] = np.zeros((n_trials, n_cells))
    coupling_coefficients[gj_conn_type] = np.zeros((n_trials, n_cells, n_cells))
    cc_conn[gj_conn_type] = []
    edge_lengths[gj_conn_type] = np.zeros((n_trials, n_cells, n_cells))
    for trial in range(n_trials):
        # load network structure
        adjacency_list = load_unique_edge_list(timestamp,
                                               gj_conn_type,
                                               trial)
        g = nx.Graph(n=n_cells)
        g.add_edges_from(adjacency_list)
        graphs[gj_conn_type].append(g)
        degrees[gj_conn_type][trial] = np.array(nx.degree(g).values())
        edge_lengths[gj_conn_type][trial] = load_edge_lengths(timestamp,
                                                              gj_conn_type,
                                                              trial)
        # load voltage responses
        responses.append(np.zeros(shape=(n_cells,n_cells), dtype=np.float))
        for stim_cell in range(n_cells):
            sim_ref = utils.cs_sim_ref(timestamp,
                                       gj_conn_type,
                                       stim_cell,
                                       trial)
            for rec_cell in range(n_cells):
                filename = '../scriptedSimulations/{0}/Golgi_network_reduced_TTX_{1}.dat'.format(sim_ref, rec_cell)
                try:
                    responses[-1][stim_cell, rec_cell] = np.loadtxt(filename)[-1]
                except IOError:
                    print('Data file not found: {0}'.format(sim_ref))
                    responses[-1][stim_cell, rec_cell] = np.NaN
        voltage_deltas = np.abs(responses[0] - responses[0].max()) #correct if hyperpolarising cell
        couplings = voltage_deltas/voltage_deltas.max(axis=1)
        couplings[np.diag_indices(n_cells)] = 0
        coupling_coefficients[gj_conn_type][trial] = couplings
        cc_conn[gj_conn_type].extend([couplings[i,j] for (i,j) in adjacency_list])
    

    cc_vs_d_ax.scatter(edge_lengths[gj_conn_type].flat,
                       coupling_coefficients[gj_conn_type].flat,
                       c=colours[gj_conn_type],
                       alpha=0.3,
                       linewidths=0)
    
    cc_ims[gj_conn_type] = cc_axes[gj_conn_type].imshow(coupling_coefficients[gj_conn_type][0], interpolation='none', cmap='coolwarm')
    cc_ims[gj_conn_type]
    cc_axes[gj_conn_type].set_title(gj_conn_type)


cc_hist_ax.hist([cc_conn['2010'], cc_conn['2012']],
                bins=10,
                color=['k', 'g'],
                label=['2010: {:.4f}$\pm${:.4f}'.format(np.mean(cc_conn['2010']),
                                                        np.std(cc_conn['2010'])),
                       '2012: {:.4f}$\pm${:.4f}'.format(np.mean(cc_conn['2012']),
                                                        np.std(cc_conn['2012']))])
cc_hist_ax.legend(loc='best')

deg_hist_ax.hist([degrees['2010'], degrees['2012']],
                bins=10,
                color=['k', 'g'],
                label=['2010: {:.1f}$\pm${:.1f}'.format(degrees['2010'].mean(), degrees['2010'].std()),
                       '2012: {:.1f}$\pm${:.1f}'.format(degrees['2012'].mean(), degrees['2012'].std())])
deg_hist_ax.legend(loc='best')

cc_fig.suptitle('Coupling coefficient')
cc_fig.subplots_adjust(bottom=0.05)
cc_cbar_ax = cc_fig.add_axes([0.15, 0.05, 0.7, 0.03])
cc_fig.colorbar(cc_ims['2012'],
                cax=cc_cbar_ax,
                orientation='horizontal')

plt.show()
