import networkx as nx
import numpy as np
import random
from matplotlib import pyplot as plt

n_nodes = 720
degree = 15

n_trials = 10
rewiring_p_range = np.logspace(-4, 0, num=12)

clustering = []
path_length = []

fig, ax = plt.subplots()
ax.set_xscale('log')
ax2 = ax.twinx()

for p in rewiring_p_range:
    print p
    c = []
    l = []
    for t in range(n_trials):
        g = nx.connected_watts_strogatz_graph(n_nodes, degree, p, tries=100)
        # for e in g.edges():
        #     g[e[0]][e[1]]['weight'] = random.gammavariate(1.53, 350)
        c.append(np.array(nx.clustering(g).values()).mean())
        l.append(nx.average_shortest_path_length(g))
    clustering.append(np.array(c).mean())
    path_length.append(np.array(l).mean())

ax.plot(rewiring_p_range, clustering)
ax2.plot(rewiring_p_range, path_length, c='k')

fig.savefig('small_world.pdf')
