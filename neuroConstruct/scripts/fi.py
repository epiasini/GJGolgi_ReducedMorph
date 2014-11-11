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

sim_config_name = 'fi_comparison'

sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
project.neuronSettings.setNoConsole()

current_amplitude_in_fa_range = range(-25000, 500, 500)


# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)

sim_refs = deque()
for amplitude_in_fa in current_amplitude_in_fa_range:
    sim_ref = 'b' + timestamp + '_' + str(int(round(amplitude_in_fa)))
    sim_refs.append(sim_ref)
    sim_path = '../simulations/' + sim_ref
    project.simulationParameters.setReference(sim_ref)
    # set current clamp amplitude
    for cell_type in ['reduced', 'vervaeke', 'solinas']:
        amplitude_in_na = amplitude_in_fa/1000000.
	stim = project.elecInputInfo.getStim('cclamp_' + cell_type)
	stim.setAmp(NumberGenerator(amplitude_in_na))
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
	print "Simulating: simulation reference " + sim_ref
	pm.doRunNeuron(sim_config)
	time.sleep(2)
	if not sim_config.getMpiConf().isRemotelyExecuted():
	    # if running locally, never have more than one sim running
	    # at the same time
	    print('Simulating on the local machine.')
	    timefile_path = sim_path + '/time.dat'
	    while not os.path.exists(timefile_path):
		time.sleep(2)

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
    time.sleep(2)

print('batch reference b' + timestamp)
System.exit(0)
