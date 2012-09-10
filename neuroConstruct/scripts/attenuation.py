import os
import random
import time
from java.lang import System
from java.io import File
from java.util import ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager

timestamp = str(time.time())
pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'attenuation'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)

project.neuronSettings.setNoConsole()

conn_name = 'relay_conn'


# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)


source_segment_index = 0
source_fraction_along = 0.5
delay = 0.
for target_segment_index in [3]:
    for target_fraction_along in [.1, .5, .9]:
	sim_ref = timestamp + '_' + str(target_segment_index) + '_' + str(target_fraction_along)
	sim_path = '../simulations/' + sim_ref
	project.simulationParameters.setReference(sim_ref)

	# delete all existing connections
	project.generatedNetworkConnections.reset()

	# connect the spike relay at the specified point on the golgi dendrite
	project.generatedNetworkConnections.addSynapticConnection(conn_name, 0, 0, source_segment_index, source_fraction_along, 0, target_segment_index, target_fraction_along, delay, None)

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
