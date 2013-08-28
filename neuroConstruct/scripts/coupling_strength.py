# Measuring the coupling strengths induced in the Golgi network by
# different gap junction connectivity rules (namely, those in
# Vervaeke2010 and Vervaeke2012).

import os
import random
import time
import csv

from java.lang import System, Long
from java.io import File
from java.util import ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.utils import NumberGenerator
#from ucl.physiol.neuroconstruct.project.cellchoice import IndividualCells

import utils

timestamp = str(time.time())
pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

n_trials = 1
sim_refs = []

for gj_conn_type in ['2010', '2012']:
    sim_config_name = 'coupling_strength_' + gj_conn_type + 'gap'
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
    project.neuronSettings.setNoConsole()
    for trial in range(n_trials):
        print('=====new trial: regenerating network=====')
        # generate
        nC_seed = Long(random.getrandbits(32))
        pm.doGenerate(sim_config_name, nC_seed)
        while pm.isGenerating():
            time.sleep(0.02)
        print('network generated')
        ##=== save network structure ===
        cell_positions_file = open(utils.cs_cell_positions_file(timestamp,
                                                                gj_conn_type,
                                                                trial),
                                   'wb')
        edge_list_file = open(utils.cs_edge_list_file(timestamp,
                                                      gj_conn_type,
                                                      trial),
                               'wb')
        cell_positions_writer = csv.writer(cell_positions_file)
        edge_list_writer = csv.writer(edge_list_file)
        # extract connectivity structure
        if gj_conn_type == '2010':
            conn_names = ['GJ2010_reduced_TTX']
        elif gj_conn_type == '2012':
            conn_names = ['GJ_reduced_TTX_GCL',
                          'GJ_reduced_TTX_ML1',
                          'GJ_reduced_TTX_ML2',
                          'GJ_reduced_TTX_ML3']
        syn_conns = [project.generatedNetworkConnections.getSynapticConnections(conn) for conn in conn_names]
        edges = []
        for conn in syn_conns:
            edges.extend([(sc.sourceEndPoint.cellNumber, sc.targetEndPoint.cellNumber) for sc in conn if sc.props[0].weight])
        # extract cell positions
        cell_positions = [(pos_record.x_pos, pos_record.y_pos, pos_record.z_pos) for pos_record in project.generatedCellPositions.getPositionRecords('Golgi_network_reduced_TTX')]
        # save to disk
        cell_positions_writer.writerows(cell_positions)
        edge_list_writer.writerows(edges)
        # close writer objects
        cell_positions_file.close()
        edge_list_file.close()

	for cell in range(45):
	    sim_ref = utils.cs_sim_ref(timestamp,
				       gj_conn_type,
				       cell,
				       trial)
	    sim_refs.append(sim_ref)
	    project.simulationParameters.setReference(sim_ref)
            ##=== simulation =====
	    # delete all existing stimuli
            #project.elecInputInfo.deleteAllStims()
	    project.generatedElecInputs.reset()

	    # set which cell to apply current clamp to
            stim = project.elecInputInfo.getStim('cclamp_network_reduced')
            #stim.setCellChooser(IndividualCells(str(cell)))
            #project.elecInputInfo.updateStim(stim)
            project.generatedElecInputs.addSingleInput(stim.getReference(),
                                                       'IClamp',
                                                       'Golgi_network_reduced_TTX',
                                                       cell,
                                                       0, 0, None)
	    # generate and compile neuron files
	    print "Generating NEURON scripts..."
	    project.neuronFileManager.setSuggestedRemoteRunTime(4)
	    simulator_seed = random.getrandbits(32)
	    project.neuronFileManager.generateTheNeuronFiles(sim_config,
							     None,
							     NeuronFileManager.RUN_HOC,
							     simulator_seed)
	    compile_process = ProcessManager(project.neuronFileManager.getMainHocFile())
	    compile_success = compile_process.compileFileWithNeuron(0,0)
	    # simulate
	    if compile_success:
		print "Submitting simulation reference " + sim_ref
		pm.doRunNeuron(sim_config)
		time.sleep(0.5) # Wait for sim to be kicked off
		if sim_config.getMpiConf().isRemotelyExecuted():
		    #utils.pull_remotes(sim_refs)
		    pass
		else:
		    # if running locally, never have more than one sim running
		    # at the same time
		    print('Simulating on the local machine.')
		    timefile_path = '../simulations/' + sim_ref + '/time.dat'
		    while not os.path.exists(timefile_path):
			time.sleep(0.5)


if sim_config.getMpiConf().isRemotelyExecuted():
    utils.wait_and_pull_remote(sim_refs, sleep_time=0.5)

print('timestamp ' + timestamp)
System.exit(0)
