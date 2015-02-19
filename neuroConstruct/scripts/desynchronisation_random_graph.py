import os
import time
import random
import networkx as nx
from math import fabs

from java.awt import Color
from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol import neuroconstruct as nc

import utils

n_trials = 32

sim_duration = 2000

simulate = True

timestamp = str(time.time())
pm = nc.project.ProjectManager(None,None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'network_desynchronisation_reduced_only_2010gap'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(sim_duration)
project.neuronSettings.setNoConsole()

remote_sim_refs = []

for trial in range(n_trials):
    # generate
    nC_seed = random.getrandbits(32)
    pm.doGenerate(sim_config_name, nC_seed)
    print('generating..')
    while pm.isGenerating():
        time.sleep(0.02)
    print('network generated')

    # extract cell positions
    cell_positions = [(pos_record.x_pos, pos_record.y_pos, pos_record.z_pos) for pos_record in project.generatedCellPositions.getAllPositionRecords()]

    sim_ref = utils.desynchronisation_random_graph(timestamp,
                                                   trial)
    project.simulationParameters.setReference(sim_ref)
    # create gap junction graph object
    gj_graph = utils.random_graph_heterogeneous_synapses(cell_positions)
    # delete all existing synaptic connections and associated information
    project.generatedNetworkConnections.reset()
    project.morphNetworkConnectionsInfo.deleteAllNetConns()
    # generate connections according to graph
    utils.nC_generate_gj_network_from_graph(project,
                                            sim_config,
                                            gj_graph,
                                            'Golgi_network_reduced',
                                            'GJGolgi_Reduced',
                                            ['GCL', 'ML1', 'ML2', 'ML3'],
                                            'Golgi_gap_2010')
    # export generated network structure to graphml for debugging
    utils.nC_network_to_graphml(project, '/home/ucbtepi/thesis/data/GoC_net_structures/graph_' + sim_ref + '.graphml')
    if simulate:
        # generate, compile and run NEURON files
        remote_sim_refs = utils.nC_generate_NEURON_and_submit(pm,
                                                              project,
                                                              sim_config,
                                                              sim_ref,
                                                              remote_sim_refs,
                                                              run_time=14)

if remote_sim_refs:
    utils.wait_and_pull_remote(remote_sim_refs, sleep_time=0.5)  

print('timestamp ' + timestamp)
System.exit(0)

