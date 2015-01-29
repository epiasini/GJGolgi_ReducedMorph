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

spatial_scale_range = [0.5, 1., 2., 128.]
n_trials = 64

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


for spatial_scale in spatial_scale_range:
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

        sim_ref = utils.variable_spatial_scale(timestamp,
                                               spatial_scale,
                                               trial)
        project.simulationParameters.setReference(sim_ref)
        # create gap junction graph object
        gj_graph = utils.spatial_graph_variable_spatial_scale(cell_positions,
                                                              spatial_scale)
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
            # generate and compile neuron files
            print "Generating NEURON scripts..."
            project.neuronFileManager.setSuggestedRemoteRunTime(14)
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
                pm.doRunNeuron(sim_config)
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

if remote_sim_refs:
    utils.wait_and_pull_remote(remote_sim_refs, sleep_time=0.5)  

print('timestamp ' + timestamp)
System.exit(0)

