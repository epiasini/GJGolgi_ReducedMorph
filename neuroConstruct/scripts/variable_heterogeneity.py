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

deg_mean_range = [9.25]#, 12.40]#range(2, 21, 1)
clustering_range = [0.52]#, 0.64]
n_trials = 1

leak_variation_fraction = 0.
sim_duration = 2000
n_cells = 45
n_cells_stimulated = 10
experimental_gjs_per_cell = 35.

timestamp = str(time.time())
pm = nc.project.ProjectManager(None,None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'variable_heterogeneity'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(sim_duration)
project.neuronSettings.setNoConsole()

remote_sim_refs = []
group_names = []

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
for i in range(n_cells):
    type_name = unicode('golgi_'+str(i))
    group_name = unicode('golgi_group_'+str(i))
    group_names.append(group_name)
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
						      'MultiDecaySyn')
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
    if i<n_cells_stimulated:
	# basolateral stimulus
	bl_delay = nc.utils.NumberGenerator()
	bl_delay.initialiseAsRandomFloatGenerator(710., 712.)#(1680., 1685.)
	bl_segment_chooser = nc.project.segmentchoice.GroupDistributedSegments('basolateral_soma', 8)
	bl_stim = nc.simulation.RandomSpikeTrainExtSettings(bl_stim_name,
							    group_name,
							    one_cell_chooser,
							    bl_segment_chooser,
							    nc.utils.NumberGenerator(0.2),
							    'MultiDecaySyn',
							    bl_delay,
							    nc.utils.NumberGenerator(10),
							    False)
	project.elecInputInfo.addStim(bl_stim)
	sim_config.addInput(bl_stim.getReference())
	# apical stimulus
	ap_delay = nc.utils.NumberGenerator()
	ap_delay.initialiseAsRandomFloatGenerator(712., 717.)#(1682., 1687.)
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
    plot = nc.project.SimPlot(type_name + '_spikes',
			      type_name + '_spikes',
			      group_name,
			      '0',
			      '0',
			      'SPIKE:-20',
			      -90,
			      50,
			      nc.project.SimPlot.SAVE_ONLY)
    project.simPlotInfo.addSimPlot(plot)
    sim_config.addPlot(type_name + '_spikes')

# generate
nC_seed = random.getrandbits(32)
pm.doGenerate(sim_config_name, nC_seed)
while pm.isGenerating():
    time.sleep(0.02)
print('network generated')

cell_positions = []
for group_name in group_names:
    positions = project.generatedCellPositions.getPositionRecords(group_name)
    assert len(positions) == 1
    cell_positions.append((positions[0].z_pos,
                           positions[0].x_pos,
                           positions[0].y_pos))

# prepare data on spatial distribution of gjs on dendritic tree
groups_with_gjs = ['GCL', 'ML1', 'ML2', 'ML3']
group_prob = (11./36, 16./36, 7./36, 2./36)
group_prob_uniform = (4./8., 2./8, 1./8., 1./8.)
group_prob = group_prob_uniform

segments_with_gjs = [golgi_reference_cell.getSegmentsInGroup(gr) for gr in groups_with_gjs]
cum_group_prob = tuple(sum(group_prob[:k]) for k in range(len(group_prob))) # cumulative

#print group_prob
#print cum_group_prob
#print segments_with_gjs
##=== set up gj network ===
for deg_mean in deg_mean_range:
    for clustering in clustering_range:
        for trial in range(n_trials):
            sim_ref = utils.variable_heterogeneity(timestamp,
                                                   deg_mean,
                                                   clustering,
                                                   trial)
            project.simulationParameters.setReference(sim_ref)
            # delete all existing synaptic connections
            project.generatedNetworkConnections.reset()
            # create gap junction graph object
            #edge_probability = float(deg_mean) / float(n_cells - 1) # p = number of edges / number of possible edges = (n_cells * deg_mean / 2) / (n_cells * (n_cells - 1) / 2)
            #gj_graph = nx.gnp_random_graph(n_cells, edge_probability)
            gj_graph = utils.clustered_poisson_graph(n_cells, deg_mean, clustering)
            #gj_graph = utils.spatial_graph_2010(cell_positions)
            # generate connections according to graph
            for i,j,data in gj_graph.edges(data=True):
                #print i, j, data
                # select random source and destination segments according to spatial distribution of gjs
                source_rand = random.random()
                source_segment_group = [k for k, p in enumerate(cum_group_prob) if p<source_rand][-1]
                source_segment = random.choice(segments_with_gjs[source_segment_group]).getSegmentId()
                dest_rand = random.random()
                dest_segment_group = [k for k, p in enumerate(cum_group_prob) if p<dest_rand][-1]
                dest_segment = random.choice(segments_with_gjs[dest_segment_group]).getSegmentId()
                #print str(i) + ',' + str(j) + ' segments ' + str(source_segment) + ',' + str(dest_segment)
                # adjust synaptic weight parameter depending on the mean of
                # the degree distribution so that in average we have the same
                # total conductance as if we had 35 gjs per cell
                #mean_gj_number = int(round(experimental_gjs_per_cell/deg_mean))
                #synaptic_weight = random.randrange(mean_gj_number-2, mean_gj_number+3, 1)
                #synaptic_weight = float(experimental_gjs_per_cell)/float(deg_mean)
                
                average_syn_weight = 516.49
                synaptic_weight = random.expovariate(1./average_syn_weight)
                #synaptic_weight = data['weight']
                conn_name = 'gj_'+str(i)+'_'+str(j)
                #synaptic_properties = nc.project.SynapticProperties('GapJuncDiscrete')
                synaptic_properties = nc.project.SynapticProperties('Golgi_gap_2010')
                synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight))
                conn_conditions = nc.project.ConnectivityConditions()
                conn_conditions.setNumConnsInitiatingCellGroup(nc.utils.NumberGenerator(0))
                project.morphNetworkConnectionsInfo.addRow(conn_name,
                                                           'golgi_group_'+str(i),
                                                           'golgi_group_'+str(j),
                                                           Vector([synaptic_properties]),
                                                           nc.project.SearchPattern.getRandomSearchPattern(),
                                                           nc.project.MaxMinLength(155., 0, 'r', 10),
                                                           conn_conditions,
                                                           Float.MAX_VALUE)
                sim_config.addNetConn(conn_name)

                connection_specific_syn_props = nc.project.ConnSpecificProps(conn_name)
                connection_specific_syn_props.weight = synaptic_weight
                project.generatedNetworkConnections.addSynapticConnection(conn_name,
                                                                          0,
                                                                          0,
                                                                          source_segment,
                                                                          random.random(),
                                                                          0,
                                                                          dest_segment,
                                                                          random.random(),
                                                                          5.,
                                                                          ArrayList([connection_specific_syn_props]))
            # export intended and resultant structure to graphml
            # file. They should coincide.
            nx.write_graphml(gj_graph, 'test_intended.graphml')
            utils.nC_network_to_graphml(project, 'test_resultant.graphml')

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

