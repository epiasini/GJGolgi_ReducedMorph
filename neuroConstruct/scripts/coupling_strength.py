# Measuring the coupling strengths induced in the Golgi network by
# different gap junction connectivity rules (namely, those in
# Vervaeke2010 and Vervaeke2012).

import os
import random
import time

from java.lang import System
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

    # generate
    pm.doGenerate(sim_config_name, 1234)
    while pm.isGenerating():
	time.sleep(0.02)
    print('network generated')

    for trial in range(n_trials):
	for cell in range(45):
	    sim_ref = utils.cs_sim_ref(timestamp,
				       gj_conn_type,
				       cell,
				       trial)
	    sim_refs.append(sim_ref)
	    project.simulationParameters.setReference(sim_ref)

	    # delete all existing stimuli
            #project.elecInputInfo.deleteAllStims()
	    project.generatedElecInputs.reset()

	    # set which cell to apply current clamp to
	    for cell_type in ['reduced']:
		stim = project.elecInputInfo.getStim('cclamp_network_' + cell_type)
		#stim.setCellChooser(IndividualCells(str(cell)))
		#project.elecInputInfo.updateStim(stim)
		project.generatedElecInputs.addSingleInput(stim.getReference(),
							   'IClamp',
							   'Golgi_network_' + cell_type + '_TTX',
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
		time.sleep(2) # Wait for sim to be kicked off
		if not sim_config.getMpiConf().isRemotelyExecuted():
		    # if running locally, never have more than one sim running
		    # at the same time
		    print('Simulating on the local machine.')
		    timefile_path = '../simulations/' + sim_ref + '/time.dat'
		    while not os.path.exists(timefile_path):
			time.sleep(5)

if sim_config.getMpiConf().isRemotelyExecuted():
    utils.wait_and_pull_remote(sim_refs, sleep_time=0.5)

print('timestamp ' + timestamp)
System.exit(0)
