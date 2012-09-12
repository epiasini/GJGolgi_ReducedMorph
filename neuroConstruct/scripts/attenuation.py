import os
import random
import time
from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.cell.utils import CellTopologyHelper

timestamp = str(time.time())
pm = ProjectManager(None, None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

sim_config_name = 'attenuation_comparison'

sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
project.neuronSettings.setNoConsole()

# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)

# pick segment ids on the detailed cell
cth = CellTopologyHelper()
cell_type = project.cellManager.getCell('GJGolgi_Reduced')
cell_type = project.cellManager.getCell('Golgi_210710_C1')
distances = dict(cth.getSegmentDistancesFromRoot(cell_type, 'all'))
dist_0 = dict((k, v) for k,v in distances.items() if 9.08 < v < 9.28)
id0 = min(dist_0.keys())
dist_1 = dict((k, v) for k,v in distances.items() if 33.78 < v < 33.98)
id1 = min(dist_1.keys())
dist_2 = dict((k, v) for k,v in distances.items() if 124.59 < v < 124.79)
id2 = min(dist_2.keys())

print id0, id1, id2

source_segment_index = 0
source_fraction_along = 0.5
delay = 0.
for conn_name in ['relay_conn', 'NetConn_relays_Golgi_Vervaeke']:
    if conn_name == 'relay_conn': target_segments = [4,5,6]
    else: target_segments = [id0, id1, id2]#[1526, 1545, 1646]#[679,1013,1196]
    for target_segment_index in target_segments:
	for target_fraction_along in [.1]:
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

print "Data reference " + timestamp + " " + str(id0) + " " + str(id1) + " " + str(id2)

System.exit(0)
