# Measuring the current-voltage relation (at soma) for detailed and
# reduced morphology cells embedded in a network with different 2010
# and 2012-style gap junctions. Simulating with different
# randomly-generated networks to test for heterogeneity.

import os
import random
import time
from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.utils import NumberGenerator

import utils

timestamp = str(time.time())
#timestamp = '1351784362.54'
pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)
timestamp_prefix = 'ir_'

stim_amplitude_range = range(-200., 220., 50.)
n_trials = 10
sim_refs = []

for gj_conn_type in ['2010', '2012']:
    sim_config_name = 'input_resistance_' + gj_conn_type + 'gap'
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
    project.neuronSettings.setNoConsole()

    for trial in range(n_trials):
	# generate
	pm.doGenerate(sim_config_name, 1234)
	while pm.isGenerating():
	    time.sleep(0.02)
	print('network generated')

	for amplitude in stim_amplitude_range:
	    sim_ref = utils.ir_sim_ref(timestamp,
				       gj_conn_type,
				       amplitude,
				       trial)
	    sim_refs.append(sim_ref)
	    sim_path = '../simulations/' + sim_ref
	    project.simulationParameters.setReference(sim_ref)
	    # set stim rate
	    amplitude_in_nA = amplitude/1000.
	    for cell_type in ['Vervaeke', 'reduced']:
		stim = project.elecInputInfo.getStim('cclamp_network_' + cell_type)
		stim.setAmp(NumberGenerator(amplitude_in_nA))
		project.elecInputInfo.updateStim(stim)
	    # generate and compile neuron files
	    print "Generating NEURON scripts..."
	    project.neuronFileManager.setSuggestedRemoteRunTime(10)
	    simulator_seed = random.getrandbits(32)
	    project.neuronFileManager.generateTheNeuronFiles(sim_config, None, NeuronFileManager.RUN_HOC,simulator_seed)
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
		    timefile_path = sim_path + '/time.dat'
		    while not os.path.exists(timefile_path):
			time.sleep(5)

if sim_config.getMpiConf().isRemotelyExecuted():
    utils.wait_and_pull_remote(sim_refs, sleep_time=5)

print('batch reference ' + timestamp_prefix + timestamp)
System.exit(0)
