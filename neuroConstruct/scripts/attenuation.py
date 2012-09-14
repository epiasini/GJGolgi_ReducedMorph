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

all_distance_limits = [[34.5, 35.5], [104.5, 105.5], [174.5, 175.5]]
reduced_seg_ids = [4,5,6]

# switch off Na channels (TTX)
reduced_cell_type = project.cellManager.getCell('GJGolgi_Reduced')
vervaeke_cell_type = project.cellManager.getCell('Golgi_210710_C1')
for chan in reduced_cell_type.getChanMechsForGroup('soma_group'):
    if chan.getName() in ['NaP_CML', 'NaR_CML', 'NaT_CML']:
	chan.setDensity(0)
	reduced_cell_type.associateGroupWithChanMech('soma_group', chan)
for chan in vervaeke_cell_type.getChanMechsForGroup('soma_group'):
    if chan.getName() in ['NaP', 'NaR', 'NaT']:
	chan.setDensity(0)
	vervaeke_cell_type.associateGroupWithChanMech('soma_group', chan)

# generate
pm.doGenerate(sim_config_name, 1234)
while pm.isGenerating():
    time.sleep(0.02)



# calculate segment-soma distances on the detailed cell
cth = CellTopologyHelper()
distances_dict = dict(cth.getSegmentDistancesFromRoot(vervaeke_cell_type, 'all'))

source_segment_index = 0
source_fraction_along = 0.5
delay = 0.

locs = []
detailed_seg_ids = []

for distance_index in [0,1,2]:
    # basic simulation setup
    sim_ref = timestamp + '_' + str(distance_index)
    sim_path = '../simulations/' + sim_ref
    project.simulationParameters.setReference(sim_ref)
    # pick segment to be stimulated on detailed cell
    dist_limits = all_distance_limits[distance_index]
    allowed_segments = dict((k,v) for k,v in distances_dict.items() if dist_limits[0] < v < dist_limits[1])
    seg_id_detailed = min(allowed_segments.keys())
    # store stimulation location for detailed cell
    dist = allowed_segments[seg_id_detailed]
    seg = vervaeke_cell_type.getSegmentWithId(seg_id_detailed)
    length = seg.getSegmentLength()
    locs.append(dist + length/2.)
    detailed_seg_ids.append(seg_id_detailed)
    # pick segment to be stimulated on reduced cell
    seg_id_reduced = reduced_seg_ids[distance_index]

    # delete all existing connections
    project.generatedNetworkConnections.reset()

    # connect the spike relay at the specified point on the golgi dendrite
    project.generatedNetworkConnections.addSynapticConnection('relay_conn',
							      0, 0, 0, 0.5, 0,
							      seg_id_reduced,
							      0.5, 0, None)
    project.generatedNetworkConnections.addSynapticConnection('NetConn_relays_Golgi_Vervaeke',
							      0, 0, 0, 0.5, 0,
							      seg_id_detailed,
							      0.5, 0, None)

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

data_string = "Data reference " + timestamp
for k,i in enumerate(detailed_seg_ids):
    data_string = data_string + " " + str(i) + " " + str(locs[k])
print data_string

System.exit(0)
