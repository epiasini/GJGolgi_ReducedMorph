import os
import random
import time
from java.lang import System
from java.io import File
from java.util import ArrayList

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

n_points = 18
detailed_distance_bounds = range(10.,190.,10.)
rec_segs_reduced = [0,4,5,6]
dendritic_group = 'apical_dend_2'

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

# pick a segment to stimulate on the detailed cell
cth = CellTopologyHelper()
distances_dict = dict(cth.getSegmentDistancesFromRoot(vervaeke_cell_type, dendritic_group))
stim_seg_detailed = random.choice([seg_id for seg_id,dist in distances_dict.items() if 200. < dist < 220.])

# pick a segment to stimulate on the reduced cell
stim_seg_reduced = 6

# map distances of all ancestors of stimulation segment
ancestor_dists = dict(cth.getDistancesFromAncestorSegments(vervaeke_cell_type, stim_seg_detailed))

source_segment_index = 0
source_fraction_along = 0.5
delay = 0.

locs = []
detailed_seg_ids = []

# basic simulation setup
sim_ref = timestamp
sim_path = '../simulations/' + sim_ref
project.simulationParameters.setReference(sim_ref)
# pick segment to be stimulated on detailed cell
allowed_segments = [[k for k,v in ancestor_dists.items() if d-2. < v < d+2.] for d in detailed_distance_bounds]
rec_segs_detailed = [stim_seg_detailed] + [random.choice(segs) for segs in allowed_segments] + [0]
# store recording distances for detailed cell
rec_dists_detailed = [0.] + [ancestor_dists[s]+vervaeke_cell_type.getSegmentWithId(s).getSegmentLength()/2.
			     for s in rec_segs_detailed[1:]]
# delete all existing connections
project.generatedNetworkConnections.reset()

# connect the spike relay at the specified point on the golgi dendrite
project.generatedNetworkConnections.addSynapticConnection('relay_conn',
							  0, 0, 0, 0.5, 0,
							  stim_seg_reduced,
							  0.5, 0, None)
project.generatedNetworkConnections.addSynapticConnection('NetConn_relays_Golgi_Vervaeke',
							  0, 0, 0, 0.5, 0,
							  stim_seg_detailed,
							  0.5, 0, None)
# delete all existing probes
project.generatedPlotSaves.reset()

# set up voltage recording at the specified points along the dendrite being stimulated
project.generatedPlotSaves.addPlotSaveDetails('Golgi_reduced_v', project.simPlotInfo.getSimPlot('Golgi_reduced_v'), ArrayList([0]), ArrayList(rec_segs_reduced), False, False)
project.generatedPlotSaves.addPlotSaveDetails('Golgi_detailed_v', project.simPlotInfo.getSimPlot('Golgi_detailed_v'), ArrayList([0]), ArrayList(rec_segs_detailed), False, False)

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
for k,i in enumerate(rec_segs_detailed):
    data_string = data_string + " " + str(i) + " " + str(rec_dists_detailed[k])
print data_string

System.exit(0)
