import os
import random
import time
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

current_amplitude_range = range(0., 440., 40.)


# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)

for amplitude in current_amplitude_range:
    sim_ref = timestamp + '_' + str(int(round(amplitude)))
    sim_path = '../simulations/' + sim_ref
    project.simulationParameters.setReference(sim_ref)
    # set current clamp amplitude
    amplitude_in_na = amplitude/1000.
    for cell_type in ['reduced', 'vervaeke', 'solinas']:
	stim = project.elecInputInfo.getStim('cclamp_' + cell_type)
	stim.setAmp(NumberGenerator(amplitude_in_na))
	project.elecInputInfo.updateStim(stim)
    # generate and compile neuron files
    print "Generating NEURON scripts..."
    simulator_seed = random.getrandbits(32)
    project.neuronFileManager.generateTheNeuronFiles(sim_config, None, NeuronFileManager.RUN_HOC,simulator_seed)
    compile_process = ProcessManager(project.neuronFileManager.getMainHocFile())
    compile_success = compile_process.compileFileWithNeuron(0,0)
    # simulate
    if compile_success:
	print "Simulating: simulation reference " + sim_ref
	pm.doRunNeuron(sim_config)
	timefile_path = sim_path + '/time.dat'
	while not os.path.exists(timefile_path):
	    time.sleep(0.1)

System.exit(0)
