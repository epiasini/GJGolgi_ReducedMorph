import os
import time
import subprocess
import random
from collections import deque
import networkx as nx
import math

def pull_remotes(sim_refs, general_sim_dir='../simulations/'):
    for sim_ref in sim_refs:
	sim_path = general_sim_dir + sim_ref
	pullsimfile_path = sim_path + '/pullsim.sh'
	print('Pulling from ' + sim_ref)
	subprocess.call([pullsimfile_path])

def wait_and_pull_remote(sim_refs, sleep_time=5, general_sim_dir='../simulations/'):
    sim_refs = deque(sim_refs)
    while sim_refs:
        # keep on going through the queue of running sims, and try to pull
	# the results from the cluster. sim_refs can a list as well as a deque.
	sim_ref = sim_refs.popleft()
	sim_path = general_sim_dir + sim_ref
	pullsimfile_path = sim_path + '/pullsim.sh'
	timefile_path = sim_path + '/time.dat'
	print('Pulling from ' + sim_ref)
	subprocess.call([pullsimfile_path])
	if not os.path.exists(timefile_path):
	    sim_refs.append(sim_ref)
	time.sleep(sleep_time)

def ir_sim_ref(timestamp, gj_conn_type, amplitude, trial):
    return 'ir' + '_' + timestamp + '_' + gj_conn_type + '_' + str(int(round(amplitude))) + '_t' + str(trial)

def ir_single_cell_sim_ref(timestamp, amplitude):
    return 'ir_sc' + '_' + timestamp + '_' + str(int(round(amplitude)))

def cs_sim_ref(timestamp, gj_conn_type, cell, trial):
    return 'cs' + '_' + timestamp + '_' + gj_conn_type + '_c' + str(cell) + '_t' + str(trial)

def cs_cell_positions_file(timestamp, gj_conn_type, trial):
    return '../scriptedSimulations/cs' + '_' + timestamp + '_' + gj_conn_type + '_t' + str(trial) + '_positions.csv'  

def cs_edge_list_file(timestamp, gj_conn_type, trial):
    return '../scriptedSimulations/cs' + '_' + timestamp + '_' + gj_conn_type + '_t' + str(trial) + '_edges.csv'  

def coupling_coefficient(r01, rl0, rl1, dv, I):
    return (rl1/(rl1+r01)) - dv*r01*(rl0+rl1+r01) / ((rl1+r01) * (I*rl0*(rl1+r01) + dv*rl0))

def variable_heterogeneity(timestamp, deg_mean, clustering, trial):
    return 'vh' + '_' + timestamp + '_m' + str(deg_mean) + '_c' + str(clustering) + '_t' + str(trial)

def poisson_cumulative(k, alpha):
    return math.exp(-alpha) * sum([alpha**i / math.factorial(i) for i in range(math.floor(k)+1)])


def poisson_sample(alpha):
    r = random.random()
    k = 0
    while True:
        if r < poisson_cumulative(k, alpha):
            break
        k += 1
    return k
    
def clustered_poisson_graph(n_nodes, mean_degree, clustering):
    is_good_triangle_sequence = False
    while not is_good_triangle_sequence:
        is_good_degree_sequence = False
        while not is_good_degree_sequence:
            degree_sequence = [poisson_sample(mean_degree) for node in range(n_nodes)]
            is_good_degree_sequence = nx.is_valid_degree_sequence(degree_sequence)
        triangle_sequence = [int(round(clustering * k * (k-1) / 2.)) for k in degree_sequence]
        is_good_triangle_sequence = (sum(triangle_sequence) % 3) == 0
        print sum(triangle_sequence)
    g = nx.random_clustered_graph(zip(degree_sequence, triangle_sequence))
    # remove parallel edges
    g = nx.Graph(g)
    # remove self loops
    g.remove_edges_from(g.selfloop_edges())
    
    return g

def distance(p, q):
    return math.sqrt(sum([(a - b)**2 for a,b in zip(p,q)]))

def connection_probability_2010(r, a=0.8223, delta=19.78, r_0=125.5):
    return a / (1 + math.exp((r-r_0)/delta))

def connection_probability_vervaeke_2010(r):
    return - 17.45 + 18.36 / (math.exp((r-267.)/39.) + 1)

def coupling_coefficient_2010(r):
    return - 2.3 + 29.7 * math.exp(-r/70.4)

def synaptic_weight_2010(r):
    cc = coupling_coefficient_2010(r)
    return 1000. * (0.576 * math.exp(cc / 12.4) + 0.000590 * math.exp(cc / 2.79) - 0.564)

def spatial_graph_2010(cell_positions,
                       connection_probability=connection_probability_vervaeke_2010,
                       synaptic_weight=synaptic_weight_2010):
    n_cells = len(cell_positions)
    edges = []
    for i, p in enumerate(cell_positions):
        for j, q in enumerate(cell_positions[i+1:]):
            d = distance(p, q)
            if random.random() < connection_probability(d):
                edges.append((i, i+1+j, {'weight': synaptic_weight(d)}))
    g = nx.Graph()
    g.add_nodes_from(range(n_cells))
    for node in g.nodes():
        g.node[node]['x'] = cell_positions[node][2]
        g.node[node]['y'] = cell_positions[node][0]
        g.node[node]['z'] = cell_positions[node][1]
    g.add_edges_from(edges)

    return g
