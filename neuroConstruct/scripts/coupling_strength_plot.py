#! /usr/bin/env python
# -*- coding: utf-8 -*-
## Coupling strength plotting script for a the Golgi network connected
## though gap junctions as defined in Vervaeke2010 or Vervaeke2012.

import sys
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from scipy.spatial import distance
from scipy.optimize import curve_fit

import utils

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
    positions = load_positions(timestamp, gj_conn_type, trial)
    distances = distance.squareform(distance.pdist(positions))
    return distances

def load_positions(timestamp, gj_conn_type, trial):
    positions = np.loadtxt(utils.cs_cell_positions_file(timestamp, gj_conn_type, trial),
                           delimiter=",")
    return positions

def exp_decay(x, a, b):
    return a*np.exp(-b * x)

def prettify_axes(ax):
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    for loc, spine in ax.spines.items():
        if loc in ['right','top']:
            spine.set_color('none') # don't draw spine

# basic controls
timestamp = sys.argv[1]
gj_conn_types = ['2010', '2012']
n_cells = 45
n_trials = 25

# load experimental data
exp_couplings = np.loadtxt('../dataSets/koen_data/golgi_pair_couplings_nonzero.csv')/100
exp_distances = np.loadtxt('../dataSets/koen_data/golgi_pair_distances_nonzero.csv')

# initialise data structures for storing results
coupling_coefficients = {} # coupling coefficients of all cell pairs
cc_conn = {} # coupling coefficients of anatomically connected cell pairs
edge_lengths = {} # pairwise distances between all cell pairs
dist_conn = {} # distances between anatomically connected cells
graphs = {} # abstract graphs representing network realisations
graphs_weighted = {} # graphs containing information on cell positions
                     # and connection strengths. These can be saved as
                     # graphml files and imported eg in Gephi for
                     # plotting.
degrees = {} # degree sequences

# initialise plotting objects
colours = {'2010':'k', '2012':'g'}
colormaps = {'2010': 'Greys', '2012': 'Greens'}
cc_axes = {}
cc_ims = {}
cc_fig, (cc_axes['2010'], cc_axes['2012']) = plt.subplots(ncols=2)
cc_vs_d_axes = {}
cc_vs_d_fig, (cc_vs_d_axes['2010'], cc_vs_d_axes['2012']) = plt.subplots(ncols=2)
deg_hist_fig, deg_hist_ax = plt.subplots()
cc_hist_fig, cc_hist_ax = plt.subplots()

for gj_conn_type in gj_conn_types:
    graphs[gj_conn_type] = []
    graphs_weighted[gj_conn_type] = []
    degrees[gj_conn_type] = np.zeros((n_trials, n_cells))
    coupling_coefficients[gj_conn_type] = np.zeros((n_trials, n_cells, n_cells))
    cc_conn[gj_conn_type] = []
    dist_conn[gj_conn_type] = []
    edge_lengths[gj_conn_type] = np.zeros((n_trials, n_cells, n_cells))
    for trial in range(n_trials):
        # load network structure
        positions = load_positions(timestamp,
                                   gj_conn_type,
                                   trial)
        adjacency_list = load_unique_edge_list(timestamp,
                                               gj_conn_type,
                                               trial)
        g = nx.Graph(n=n_cells)
        g.add_edges_from(adjacency_list)
        graphs[gj_conn_type].append(g)
        graphs_weighted[gj_conn_type].append(g)
        for node in range(n_cells):
            g.node[node]['x'] = positions[node][0]
            g.node[node]['y'] = positions[node][2]
        degrees[gj_conn_type][trial] = np.array(nx.degree(g).values())
        edge_lengths[gj_conn_type][trial] = load_edge_lengths(timestamp,
                                                              gj_conn_type,
                                                              trial)
        # load voltage responses
        responses = np.zeros(shape=(n_cells,n_cells), dtype=np.float)
        for stim_cell in range(n_cells):
            sim_ref = utils.cs_sim_ref(timestamp,
                                       gj_conn_type,
                                       stim_cell,
                                       trial)
            for rec_cell in range(n_cells):
                filename = '../scriptedSimulations/{0}/Golgi_network_reduced_TTX_{1}.dat'.format(sim_ref, rec_cell)
                try:
                    responses[stim_cell, rec_cell] = np.loadtxt(filename)[-1]
                except IOError:
                    print('Data file not found: {0}'.format(sim_ref))
                    responses[stim_cell, rec_cell] = np.NaN
        ##=== trial-specific analysis ===
        voltage_deltas = np.abs(responses - responses.max()) #correct if hyperpolarising cell and if there is at least a pair of cells with near-zero coupling
        couplings = voltage_deltas/voltage_deltas.max(axis=1)
        couplings[np.diag_indices(n_cells)] = 0
        coupling_coefficients[gj_conn_type][trial] = couplings
        cc_conn[gj_conn_type].extend([couplings[i,j] for (i,j) in adjacency_list])
        dist_conn[gj_conn_type].extend([edge_lengths[gj_conn_type][trial][i,j] for (i,j) in adjacency_list])
        for (i,j) in adjacency_list:
            graphs_weighted[gj_conn_type][trial][i][j]['weight'] = couplings[i,j]

    ##=== connection type-specific analysis ===
    # only include cell pairs whose distance is less than 160um when
    # estimating the "cc vs distance" relationship
    # cell_idxs = edge_lengths[gj_conn_type].flat > 0
    # cc_vs_d_axes[gj_conn_type].scatter(edge_lengths[gj_conn_type].flat[cell_idxs],
    #                                    coupling_coefficients[gj_conn_type].flat[cell_idxs],
    #                                    c=colours[gj_conn_type],
    #                                    alpha=0.2,
    #                                    linewidths=0)
    cc_vs_d_axes[gj_conn_type].scatter(dist_conn[gj_conn_type],
                                       cc_conn[gj_conn_type],
                                       c=colours[gj_conn_type],
                                       alpha=0.2,
                                       linewidths=0)
    a, b = curve_fit(exp_decay,
                     np.array(dist_conn[gj_conn_type]),
                     np.array(cc_conn[gj_conn_type]),
                     p0=[0.5, 1./60])[0]
    fit_x_model = np.arange(10, 160, 0.1)
    cc_vs_d_axes[gj_conn_type].plot(fit_x_model,
                                    exp_decay(fit_x_model, a, b),
                                    c='y',
                                    alpha=1,
                                    linewidth=1.5)
    fit_x_exp = np.arange(10, 160, 0.1)
    cc_vs_d_axes[gj_conn_type].scatter(exp_distances,
                                       exp_couplings,
                                       c='r',
                                       alpha=0.5,
                                       linewidths=0)
    cc_vs_d_axes[gj_conn_type].plot(fit_x_exp,
                                    0.01*(-2.3+29.7*np.exp(-fit_x_exp/70.4)),
                                    c='r',
                                    linewidth=1.5)
    prettify_axes(cc_vs_d_axes[gj_conn_type])
    cc_vs_d_axes[gj_conn_type].set_xlabel(r'distance ($\mu m$)')
    cc_vs_d_axes[gj_conn_type].set_ylabel('coupling coefficient')
    # display the coupling coefficient dissimilarity matrix for a
    # sample instantiation of this connectivity model
    cc_ims[gj_conn_type] = cc_axes[gj_conn_type].imshow(coupling_coefficients[gj_conn_type][0],
                                                        interpolation='none',
                                                        cmap='bone',
                                                        vmax=0.5)
    cc_ims[gj_conn_type]
    cc_axes[gj_conn_type].set_title(gj_conn_type)
    # save a sample instantiation of this connectivity model as a
    # weighted graph where edges are present only between directly
    # connected cells and edge weights corresopnd to coupling
    # coefficients.
    nx.write_graphml(graphs_weighted[gj_conn_type][0],
                     '{}_weighted.graphml'.format(gj_conn_type))

##=== global analysis ===



##=== plotting ===
cc_hist_ax.hist([cc_conn['2010'], cc_conn['2012']],
                bins=15,
                normed=True,
                histtype='bar',
                color=['k', 'g'],
                alpha=0.6,
                label=['2010: {:.4f}$\pm${:.4f}'.format(np.mean(cc_conn['2010']),
                                                        np.std(cc_conn['2010'])),
                       '2012: {:.4f}$\pm${:.4f}'.format(np.mean(cc_conn['2012']),
                                                        np.std(cc_conn['2012']))])
cc_hist_ax.set_xlabel('coupling coefficient')
cc_hist_fig.suptitle('Coupling coefficient for directly connected pairs')
cc_hist_ax.legend(loc='best')
prettify_axes(cc_hist_ax)

min_degree = min(degrees['2010'].min(), degrees['2012'].min())
max_degree = max(degrees['2010'].max(), degrees['2012'].max())
deg_hist_ax.hist([degrees['2010'], degrees['2012']],
                bins=max_degree-min_degree+1,
                range=(min_degree-0.5, max_degree+0.5),
                histtype='bar',
                normed=True,
                color=['k', 'g'],
                alpha=0.6,
                label=['2010: {:.1f}$\pm${:.1f}'.format(degrees['2010'].mean(), degrees['2010'].std()),
                       '2012: {:.1f}$\pm${:.1f}'.format(degrees['2012'].mean(), degrees['2012'].std())])
deg_hist_ax.set_xlabel('degree')
deg_hist_fig.suptitle('Degree distribution')
deg_hist_ax.legend(loc='best')
prettify_axes(deg_hist_ax)

cc_fig.suptitle('Coupling coefficient - single network realisation')
cc_fig.subplots_adjust(bottom=0.05)
cc_cbar_ax = cc_fig.add_axes([0.15, 0.05, 0.7, 0.03])
cc_fig.colorbar(cc_ims['2012'],
                cax=cc_cbar_ax,
                orientation='horizontal')

plt.show()
