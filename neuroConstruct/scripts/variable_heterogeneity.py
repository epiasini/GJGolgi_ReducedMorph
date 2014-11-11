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

deg_mean_range = [7]
deg_sigma_cv_range = [.35]
leak_variation_fraction = 0.2
sim_duration = 4000

timestamp = str(time.time())
pm = nc.project.ProjectManager(None,None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'variable_heterogeneity'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(sim_duration)
project.neuronSettings.setNoConsole()

## ==== introduce single-cell-level heterogeneity in somatic leak conductance ===
# boilerplate stuff
golgi_reference_cell = project.cellManager.getCell('GJGolgi_Reduced')
#golgi_reference_cell = project.cellManager.getCell('Golgi_210710_C1')
region_name = 'Regions_1'
leak_cond_name = 'LeakCond'
colour = Color(255, 51, 51)
one_cell_chooser = nc.project.cellchoice.FixedNumberCells(1)
adapter = nc.project.packing.RandomCellPackingAdapter()
adapter.setParameter(nc.project.packing.RandomCellPackingAdapter.CELL_NUMBER_POLICY, 1)
# subtract maximum variation in leak conductance from base value
base_golgi_leak = [c.getDensity() for c in golgi_reference_cell.getChanMechsForGroup('all') if c.getName()==leak_cond_name][0]
max_leak_cond_delta = base_golgi_leak*leak_variation_fraction
new_golgi_leak = base_golgi_leak - max_leak_cond_delta
for chan in golgi_reference_cell.getChanMechsForGroup('all'):
    if chan.getName() == leak_cond_name:
	chan.setDensity(new_golgi_leak)
	golgi_reference_cell.associateGroupWithChanMech('all', chan)
for i in range(45):
    type_name = unicode('golgi_'+str(i))
    group_name = unicode('golgi_group_'+str(i))
    bl_noise_name = unicode('golgi_noise_'+str(i)+'_bl')
    ap_noise_name = unicode('golgi_noise_'+str(i)+'_ap')
    bl_stim_name = unicode('golgi_stim_'+str(i)+'_bl')
    ap_stim_name = unicode('golgi_stim_'+str(i)+'_ap')
    # create and add new cell type
    new_cell_type = golgi_reference_cell.clone()
    new_cell_type.setInstanceName(type_name)
    project.cellManager.addCellType(new_cell_type)
    # create and add new cell group
    project.cellGroupsInfo.addCellGroup(group_name,
                                        type_name,
                                        region_name,
                                        colour,
                                        adapter,
                                        i)
    sim_config.addCellGroup(group_name)
    # set cell-specific value for somatic leak conductance
    for chan in new_cell_type.getChanMechsForGroup('all'):
	if chan.getName() == 'VariableLeakConductance':
	    chan.setDensity(random.uniform(0.,
					   2.*max_leak_cond_delta))
	    new_cell_type.associateGroupWithChanMech('all', chan)
    ## =create and add new stimuli=
    # basolateral background
    bl_noise_segchooser = nc.project.segmentchoice.GroupDistributedSegments('basolateral_soma', 20)
    bl_noise = nc.simulation.RandomSpikeTrainSettings(bl_noise_name,
						      group_name,
						      one_cell_chooser,
						      bl_noise_segchooser,
						      nc.utils.NumberGenerator(0.002),
						      'Golgi_AMPA_mf')
    project.elecInputInfo.addStim(bl_noise)
    sim_config.addInput(bl_noise.getReference())
    # apical background
    ap_noise_segchooser = nc.project.segmentchoice.GroupDistributedSegments('parallel_fibres', 100)
    ap_noise = nc.simulation.RandomSpikeTrainSettings(ap_noise_name,
						      group_name,
						      one_cell_chooser,
						      ap_noise_segchooser,
						      nc.utils.NumberGenerator(0.0005),
						      'ApicalSyn')
    project.elecInputInfo.addStim(ap_noise)
    sim_config.addInput(ap_noise.getReference())
    if i<5:
	# basolateral stimulus
	bl_delay = nc.utils.NumberGenerator()
	bl_delay.initialiseAsRandomFloatGenerator(1680., 1685.)
	bl_segment_chooser = nc.project.segmentchoice.GroupDistributedSegments('basolateral_soma', 8)
	bl_stim = nc.simulation.RandomSpikeTrainExtSettings(bl_stim_name,
							    group_name,
							    one_cell_chooser,
							    bl_segment_chooser,
							    nc.utils.NumberGenerator(0.2),
							    'Golgi_AMPA_mf',
							    bl_delay,
							    nc.utils.NumberGenerator(10),
							    False)
	project.elecInputInfo.addStim(bl_stim)
	sim_config.addInput(bl_stim.getReference())
	# apical stimulus
	ap_delay = nc.utils.NumberGenerator()
	ap_delay.initialiseAsRandomFloatGenerator(1682., 1687.)
	ap_segment_chooser = nc.project.segmentchoice.GroupDistributedSegments('parallel_fibres', 50)
	ap_stim = nc.simulation.RandomSpikeTrainExtSettings(ap_stim_name,
							    group_name,
							    one_cell_chooser,
							    ap_segment_chooser,
							    nc.utils.NumberGenerator(0.35),
							    'ApicalSyn',
							    ap_delay,
							    nc.utils.NumberGenerator(15),
							    False)

	project.elecInputInfo.addStim(ap_stim)
	sim_config.addInput(ap_stim.getReference())
    # create and add new plot/save objects
    plot = nc.project.SimPlot(type_name + '_v',
			      type_name + '_v',
			      group_name,
			      '0',
			      '0',
			      nc.project.SimPlot.VOLTAGE,
			      -90,
			      50,
			      nc.project.SimPlot.SAVE_ONLY)
    project.simPlotInfo.addSimPlot(plot)
    sim_config.addPlot(type_name + '_v')

# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)
print('network generated')

# prepare data on spatial distribution of gjs on dendritic tree
groups_with_gjs = ['GCL', 'ML1', 'ML2', 'ML3']
group_prob = (11./36, 16./36, 7./36, 2./36)
segments_with_gjs = [golgi_reference_cell.getSegmentsInGroup(gr) for gr in groups_with_gjs]
cum_group_prob = tuple(sum(group_prob[:k]) for k in range(len(group_prob))) # cumulative

##=== set up gj network ===
for deg_mean in deg_mean_range:
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
	    # select random source and destination segments according to spatial distribution of gjs
	    source_rand = random.random()
	    source_segment_group = [k for k, p in enumerate(cum_group_prob) if source_rand>=p][-1]
	    source_segment = random.choice(segments_with_gjs[source_segment_group]).getSegmentId()
	    dest_rand = random.random()
	    dest_segment_group = [k for k, p in enumerate(cum_group_prob) if dest_rand>=p][-1]
	    dest_segment = random.choice(segments_with_gjs[dest_segment_group]).getSegmentId()
	    print str(i) + ',' + str(j) + ' segments ' + str(source_segment) + ',' + str(dest_segment)
	    # adjust synaptic weight parameter depending on the mean
	    # of the degree distribution so that in average we have 35
	    # gjs per cell, but allowing for some variability on the
	    # number of gjs per connected pair.
	    mean_gj_number = int(round(35./deg_mean))
	    synaptic_weight = random.randrange(mean_gj_number-2, mean_gj_number+3, 1)
	    conn_name = 'gj_'+str(i)+'_'+str(j)
	    synaptic_properties = nc.project.SynapticProperties('GapJuncDiscrete')
	    synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight))
	    synaptic_properties_list = Vector([synaptic_properties])
            conn_conditions = nc.project.ConnectivityConditions()
            conn_conditions.setNumConnsInitiatingCellGroup(nc.utils.NumberGenerator(0))
	    project.morphNetworkConnectionsInfo.addRow(conn_name,
						       'golgi_group_'+str(i),
						       'golgi_group_'+str(j),
						       synaptic_properties_list,
						       nc.project.SearchPattern.getRandomSearchPattern(),
						       nc.project.MaxMinLength(Float.MAX_VALUE, 0, 'r', 100),
						       conn_conditions,
						       Float.MAX_VALUE)
	    sim_config.addNetConn(conn_name)
	    project.generatedNetworkConnections.addSynapticConnection(conn_name,
								      0,
								      0,
								      source_segment,
								      0.5,
								      0,
								      dest_segment,
								      0.5,
								      0,
								      None)
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
