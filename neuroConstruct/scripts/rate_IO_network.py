# Sets up a series of simulations to estimate the I/O rate-coded
# relationship for various Golgi cell network models. Of course, what
# this script does depends on the state of the neuroConstruct project.
import os
import random
import time
import subprocess
from collections import deque
from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.cell.utils import CellTopologyHelper
from ucl.physiol.neuroconstruct.utils import NumberGenerator

timestamp = str(time.time())
pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

timestamp_prefix = 'net'
stim_source = 'mf'

sim_config_name = 'network_IO_' + stim_source
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
project.neuronSettings.setNoConsole()

stim_rate_range = range(1., 440., 40.) # neuron doesn't like zero frequency


# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)
print('network generated')

sim_refs = deque()
for rate in stim_rate_range:
    sim_ref = timestamp_prefix + timestamp + '_' + str(int(round(rate)))
    sim_refs.append(sim_ref)
    sim_path = '../simulations/' + sim_ref
    project.simulationParameters.setReference(sim_ref)
    # set stim rate
    rate_in_kHz = rate/1000.
    for cell_type in ['Vervaeke', 'reduced']:
	stim = project.elecInputInfo.getStim('network_' + cell_type + '_' +
					     stim_source + '_AMPAonly')
	stim.setRate(NumberGenerator(rate_in_kHz))
	project.elecInputInfo.updateStim(stim)
    # generate and compile neuron files
    print "Generating NEURON scripts..."
    project.neuronFileManager.setSuggestedRemoteRunTime(20)
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

while sim_refs and sim_config.getMpiConf().isRemotelyExecuted():
    # keep on going through the queue of running sims, and try to pull
    # the results from the cluster
    sim_ref = sim_refs.popleft()
    sim_path = '../simulations/' + sim_ref
    pullsimfile_path = sim_path + '/pullsim.sh'
    timefile_path = sim_path + '/time.dat'
    print('Pulling from ' + sim_ref)
    subprocess.call([pullsimfile_path])
    if not os.path.exists(timefile_path):
	sim_refs.append(sim_ref)
    time.sleep(5)

print('batch reference net' + timestamp)
System.exit(0)
