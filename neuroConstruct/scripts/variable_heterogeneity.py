import os
import time
import random
import networkx as nx
from math import fabs

from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol import neuroconstruct as nc

import utils

deg_mean_range = [35]
deg_sigma_cv_range = [1./3]
syn_strength_noise = 0.6

timestamp = str(time.time())
pm = nc.project.ProjectManager(None,None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'variable_heterogeneity'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
project.neuronSettings.setNoConsole()
# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)
print('network generated')

for deg_mean in deg_mean_range:
    # adjust synaptic strenght to keep average coupling conductance constant
    deg_sigma_range = [x*deg_mean for x in deg_sigma_cv_range]

    for deg_sigma in deg_sigma_range:
	sim_ref = utils.variable_heterogeneity(timestamp,
					       deg_mean,
					       deg_sigma)
	project.simulationParameters.setReference(sim_ref)
	# delete all existing connections
	project.generatedNetworkConnections.reset()
	project.morphNetworkConnectionsInfo.deleteAllNetConns()
        sim_config.setNetConns(ArrayList())
	# generate a degree sequence from the appropriate distribution
	deg_seq = []
	is_good_sequence = False
	while not is_good_sequence:
	    deg_seq = [int(fabs(round(random.gauss(mu=deg_mean,sigma=deg_sigma)))) for each in range(45)]
	    is_good_sequence = nx.is_valid_degree_sequence(deg_seq)
	# create gap junction graph object
	gj_graph = nx.configuration_model(deg_seq)
	# remove self edges
	gj_graph.remove_edges_from(gj_graph.selfloop_edges())
	# generate connections according to graph
	for i,j in gj_graph.edges():
	    conn_name = 'gj_'+str(i)+'_'+str(j)
	    synaptic_weight = 35./deg_mean
	    synaptic_properties = nc.project.SynapticProperties('GapJuncDiscrete')
	    synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight*(1+(random.random()-0.5)*syn_strength_noise)))
	    synaptic_properties_list = Vector([synaptic_properties])
            conn_conditions = nc.project.ConnectivityConditions()
            conn_conditions.setNumConnsInitiatingCellGroup(nc.utils.NumberGenerator(0))
	    project.morphNetworkConnectionsInfo.addRow(conn_name, 'Golgi_network_reduced', 'Golgi_network_reduced', synaptic_properties_list, nc.project.SearchPattern.getRandomSearchPattern(), nc.project.MaxMinLength(Float.MAX_VALUE, 0, 'r', 100), conn_conditions, Float.MAX_VALUE)
	    sim_config.addNetConn(conn_name)
	    project.generatedNetworkConnections.addSynapticConnection(conn_name, i, j)
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
	    time.sleep(0.5) # Wait for sim to be kicked off
	    print('Simulating on the local machine.')
	    timefile_path = '../simulations/' + sim_ref + '/time.dat'
	    while not os.path.exists(timefile_path):
		time.sleep(0.5)

print('timestamp ' + timestamp)
System.exit(0)

