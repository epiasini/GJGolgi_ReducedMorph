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

def cs_sim_ref(gj_conn_type, variance_scaling, cell, trial):
    return 'cs' + '_' + gj_conn_type + '_c' + str(cell) + '_t' + str(trial)

def coupling_coefficient(r01, rl0, rl1, dv, I):
    return (rl1/(rl1+r01)) - dv*r01*(rl0+rl1+r01) / ((rl1+r01) * (I*rl0*(rl1+r01) + dv*rl0))

def variable_heterogeneity(timestamp, mean_scaling, variance_scaling, trial):
    return 'vh' + '_' + timestamp + '_ms' + str(mean_scaling) + '_vs' + str(variance_scaling) + '_t' + str(trial)

def variable_spatial_scale(timestamp, spatial_scale, trial):
    return 'vs_' + timestamp + '_ss' + str(spatial_scale) + '_t' + str(trial)

def desynchronisation_small_world(timestamp, rewiring_p, trial):
    return 'sw_' + timestamp + '_rw' + str(rewiring_p) + '_t' + str(trial)

def desynchronisation_random_graph(timestamp, trial):
    return 'rg_' + timestamp + '_t' + str(trial)


def distance(p, q):
    return math.sqrt(sum([(a - b)**2 for a,b in zip(p,q)]))

def connection_probability_2010(r, a=0.8223, delta=19.78, r_0=125.5):
    return a / (1 + math.exp((r-r_0)/delta))

def connection_probability_vervaeke_2010(r):
    return - 17.45 + 18.36 / (math.exp((r-267.)/39.) + 1)

def coupling_coefficient_vervaeke_2010(r):
    return - 2.3 + 29.7 * math.exp(-r/70.4)

def synaptic_weight_vervaeke_2010(r):
    cc = coupling_coefficient_vervaeke_2010(r)
    return 1000. * (0.576 * math.exp(cc / 12.4) + 0.000590 * math.exp(cc / 2.79) - 0.564)

def spatial_graph_2010(cell_positions,
                       connection_probability=connection_probability_vervaeke_2010,
                       synaptic_weight=synaptic_weight_vervaeke_2010,
                       state=None):
    if state:
        random.setstate(state)
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
        g.node[node]['x'] = cell_positions[node][0]
        g.node[node]['y'] = cell_positions[node][1]
        g.node[node]['z'] = cell_positions[node][2]
    g.add_edges_from(edges)
    return g

def spatial_graph_variable_spatial_scale(cell_positions,
                                         spatial_scale=1.,
                                         connection_probability=connection_probability_vervaeke_2010,
                                         synaptic_weight=synaptic_weight_vervaeke_2010):
    state = random.getstate()
    g_2010 = spatial_graph_2010(cell_positions)
    weights_2010 = [e[2]['weight'] for e in g_2010.edges(data=True)]
    total_weight_2010 = sum(weights_2010)
    # reset RNG to make sure we will rescale strengths fairly
    random.setstate(state)

    # generate spatial network with 2010 rules but scaling all distances
    n_cells = len(cell_positions)
    edges = []
    for i, p in enumerate(cell_positions):
        for j, q in enumerate(cell_positions[i+1:]):
            d = distance(p, q) / spatial_scale
            if random.random() < connection_probability(d):
                edges.append((i, i+1+j, {'weight': synaptic_weight(d)}))

    # rescale weights to keep the same total value across the network
    weights = [e[2]['weight'] for e in edges]
    total_weight = sum(weights)
    for e in edges:
        e[2]['weight'] *= total_weight_2010 / total_weight

    # create graph object
    g = nx.Graph()
    g.add_nodes_from(range(n_cells))
    for node in g.nodes():
        g.node[node]['x'] = cell_positions[node][0]
        g.node[node]['y'] = cell_positions[node][1]
        g.node[node]['z'] = cell_positions[node][2]
    g.add_edges_from(edges)
    return g

def spatial_graph_arbitrary_variance(cell_positions, mean_scaling=1., variance_scaling=1.):
    # create graph with same spatial structure as 2010 model
    g = spatial_graph_2010(cell_positions)
    # define k and theta values for distribution of synaptic weights
    # fitted on many realisations of the 2010 model
    fitted_k = 1.53
    fitted_theta = 350.
    # scale k and theta according to scaling parameter
    k = fitted_k * mean_scaling / variance_scaling
    theta = fitted_theta * variance_scaling
    # set synaptic weights according to gamma(k,theta) 
    for e in g.edges():
        g[e[0]][e[1]]['weight'] = random.gammavariate(k, theta)
    return g

def spatial_graph_shuffled_weights(cell_positions):
    g = spatial_graph_2010(cell_positions)
    weights = [e[2]['weight'] for e in g.edges(data=True)]
    # shuffle edge weights
    random.shuffle(weights)
    for k, e in enumerate(g.edges()):
        g[e[0]][e[1]]['weight'] = weights[k]
    return g

def random_graph_heterogeneous_synapses(cell_positions):
    h = spatial_graph_2010(cell_positions)
    weights = [e[2]['weight'] for e in h.edges(data=True)]
    random.shuffle(weights)
    g = nx.gnm_random_graph(h.order(), h.size())
    for k, e in enumerate(g.edges()):
        g[e[0]][e[1]]['weight'] = weights[k]
    return g

def small_world_graph(n_cells, degree, rewiring_p, tries=100):
    g = nx.connected_watts_strogatz_graph(n_cells, degree, rewiring_p, tries)
    for e in g.edges():
        g[e[0]][e[1]]['weight'] = random.gammavariate(1.53, 350)
    return g
    
def nC_network_to_graphml(project, graphml_file_path='test.graphml'):
    # extract cell position records
    cell_positions = [(pos_record.x_pos, pos_record.y_pos, pos_record.z_pos) for pos_record in project.generatedCellPositions.getAllPositionRecords()]

    # create graph object
    graph = nx.Graph()
    graph.add_nodes_from(range(len(cell_positions)))

    # add node properties for positions. Permute x,y,z to get a nicer
    # default visualisation in Gephi when opening the resulting
    # graphml file.
    for k, position in enumerate(cell_positions):
        graph.node[k]['x'] = position[2]
        graph.node[k]['y'] = position[0]
        graph.node[k]['z'] = position[1]

    # add edges
    for conn_name in project.generatedNetworkConnections.getNamesNonEmptyNetConns():
        conns = project.generatedNetworkConnections.getSynapticConnections(conn_name)
        assert len(conns) == 1
        conn = conns[0]
        source = conn.sourceEndPoint.cellNumber
        target = conn.targetEndPoint.cellNumber
        graph.add_edge(source, target, weight=conn.props[0].weight)

    # save to disk
    nx.write_graphml(graph, graphml_file_path)
    
    return graph

def nC_generate_gj_network_from_graph(project,
                                      sim_config,
                                      graph,
                                      cell_group_name,
                                      cell_model_name,
                                      segment_group_names,
                                      synapse_model_name):
    from java.lang import Float
    from java.util import Vector, ArrayList
    from ucl.physiol import neuroconstruct as nc

    cell_model = project.cellManager.getCell(cell_model_name)
    dendritic_segments = []
    for gr in segment_group_names:
        dendritic_segments.extend(cell_model.getSegmentsInGroup(gr))

    for source_cell, target_cell, data in graph.edges(data=True):
        synaptic_weight = data['weight']
        # select random source and destination segments for
        # synaptic connection
        source_segment = random.choice(dendritic_segments).getSegmentId()
        target_segment = random.choice(dendritic_segments).getSegmentId()

        conn_name = 'gj_'+str(source_cell)+'_'+str(target_cell)

        synaptic_properties = nc.project.SynapticProperties(synapse_model_name)
        synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight))
        conn_conditions = nc.project.ConnectivityConditions()
        conn_conditions.setNumConnsInitiatingCellGroup(nc.utils.NumberGenerator(0))
        project.morphNetworkConnectionsInfo.addRow(conn_name,
                                                   cell_group_name,
                                                   cell_group_name,
                                                   Vector([synaptic_properties]),
                                                   nc.project.SearchPattern.getRandomSearchPattern(),
                                                   nc.project.MaxMinLength(Float.MAX_VALUE, 0, 'r', 10),
                                                   conn_conditions,
                                                   Float.MAX_VALUE)
        sim_config.addNetConn(conn_name)

        connection_specific_syn_props = nc.project.ConnSpecificProps(conn_name)
        connection_specific_syn_props.weight = synaptic_weight
        project.generatedNetworkConnections.addSynapticConnection(conn_name,
                                                                  0,
                                                                  source_cell,
                                                                  source_segment,
                                                                  random.random(),
                                                                  target_cell,
                                                                  target_segment,
                                                                  random.random(),
                                                                  5.,
                                                                  ArrayList([connection_specific_syn_props]))
        

def nC_generate_NEURON_and_submit(project_manager,
                                  project,
                                  sim_config,
                                  sim_ref,
                                  remote_sim_refs,
                                  run_time=60):
    from ucl.physiol import neuroconstruct as nc

    print "Generating NEURON scripts..."
    project.neuronFileManager.setSuggestedRemoteRunTime(run_time)
    simulator_seed = random.getrandbits(32)
    project.neuronFileManager.generateTheNeuronFiles(sim_config,
                                                     None,
                                                     nc.neuron.NeuronFileManager.RUN_HOC,
                                                     simulator_seed)
    compile_process = nc.nmodleditor.processes.ProcessManager(project.neuronFileManager.getMainHocFile())
    compile_success = compile_process.compileFileWithNeuron(0,0)
    # simulate
    if compile_success:
        print "Submitting simulation reference " + sim_ref
        project_manager.doRunNeuron(sim_config)
        time.sleep(2.5) # Wait for sim to be kicked off
        if sim_config.getMpiConf().isRemotelyExecuted():
            remote_sim_refs.append(sim_ref)
        else:
            # if running locally, never have more than one sim running
            # at the same time
            print('Simulating on the local machine.')
            timefile_path = '../simulations/' + sim_ref + '/time.dat'
            while not os.path.exists(timefile_path):
                time.sleep(0.5)

    return remote_sim_refs
