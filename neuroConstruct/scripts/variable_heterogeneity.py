import os
import time
import random
import networkx as nx
from math import fabs

from java.lang import System
from java.io import File
from java.util import Vector

from ucl.physiol import neuroconstruct as nc

import utils

deg_mean_range = [10,35]
deg_sigma_cv_range = [1./20, 1./4, 1./2]

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
    synaptic_weight = 35./deg_mean
    synaptic_properties = nc.project.SynapticProperties('GapJuncDiscrete')
    synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight))
    synaptic_properties_list = Vector([synaptic_properties])
    project.morphNetworkConnectionsInfo.setSynapseList('variable_heterogeneity_gj',
						   synaptic_properties_list)
    deg_sigma_range = [x*deg_mean for x in deg_sigma_cv_range]

    for deg_sigma in deg_sigma_range:
	sim_ref = utils.variable_heterogeneity(timestamp,
					       deg_mean,
					       deg_sigma)
	project.simulationParameters.setReference(sim_ref)
	# delete all existing connections
	project.generatedNetworkConnections.reset()
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
	    project.generatedNetworkConnections.addSynapticConnection('variable_heterogeneity_gj',
								      i,
								      j)
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

