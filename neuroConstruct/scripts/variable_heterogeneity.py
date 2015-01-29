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

variance_scaling_range = [0.125, 8.]
n_trials = 1

sim_duration = 700

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


cell_model = project.cellManager.getCell('GJGolgi_Reduced')


# prepare data on spatial distribution of gjs on dendritic tree
dendritic_segments = []
for gr in ['GCL', 'ML1', 'ML2', 'ML3']:
    dendritic_segments.extend(cell_model.getSegmentsInGroup(gr))

for variance_scaling in variance_scaling_range:
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
        
        sim_ref = utils.variable_heterogeneity(timestamp,
                                               variance_scaling,
                                               trial)
        project.simulationParameters.setReference(sim_ref)
        # delete all existing synaptic connections
        project.generatedNetworkConnections.reset()
        # create gap junction graph object
        gj_graph = utils.spatial_graph_arbitrary_variance(cell_positions,
                                                          variance_scaling)
        # generate connections according to graph
        for source_cell, target_cell, data in gj_graph.edges(data=True):
            # select random source and destination segments for
            # synaptic connection
            source_segment = random.choice(dendritic_segments).getSegmentId()
            target_segment = random.choice(dendritic_segments).getSegmentId()

            conn_name = 'GJ2010_reduced'
            connection_specific_syn_props = nc.project.ConnSpecificProps(conn_name)
            connection_specific_syn_props.weight = data['weight']
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
        # export generated network structure to graphml for debugging
        utils.nC_network_to_graphml(project,
                                    'GJ2010_reduced',
                                    '/home/ucbtepi/thesis/data/GoC_net_structures/graph_randomised_vs' + str(variance_scaling) + '_trial' + str(trial) +'.graphml')
        if simulate:
            # generate and compile neuron files
            print "Generating NEURON scripts..."
            project.neuronFileManager.setSuggestedRemoteRunTime(40)
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

