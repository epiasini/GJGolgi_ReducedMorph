import os
import random
import time
from java.lang import System
from java.io import File
from java.util import ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager


pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'attenuation'
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)

conn_name = 'relay_conn'
source_segment_index = 0
source_fraction_along = 0.5
target_segment_index = 2
target_fraction_along = 0.5
delay = 0.


# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)

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
    print "Simulating"
    pm.doRunNeuron(sim_config)

System.exit(0)
