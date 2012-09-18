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
detailed_distance_step = 180/n_points
detailed_distance_bounds = [[9.5, 10.5], [19.5, 20.5], [29.5,30.5],
			    [39.5, 40.5], [49.5, 50.5], [59.5, 60.5],
			    [69.5, 70.5], [79.5, 80.5], [89.5, 90.5],
			    [99.5, 100.5], [109.5, 110.5], [119.5, 120.5],
			    [129.5, 130.5], [139.5, 140.5], [149.5, 150.5],
			    [159.5, 160.5], [169.5, 170.5], [179.5, 180.5]]
reduced_seg_ids = [4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6]

# switch off Na channels (TTX)
reduced_cell_type = project.cellManager.getCell('GJGolgi_Reduced')
vervaeke_cell_type = project.cellManager.getCell('Golgi_240710_C0')
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

for distance_index in range(n_points):
    # basic simulation setup
    sim_ref = timestamp + '_' + str(distance_index)
    sim_path = '../simulations/' + sim_ref
    project.simulationParameters.setReference(sim_ref)
    # pick segment to be stimulated on detailed cell
    dist_limits = detailed_distance_bounds[distance_index]
    allowed_segments = dict((k,v) for k,v in distances_dict.items() if dist_limits[0] < v < dist_limits[1])
    seg_id_detailed = random.choice(allowed_segments.keys())
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
    # delete all existing probes
    project.generatedPlotSaves.reset()

    # set up voltage recording at the soma
    project.generatedPlotSaves.addPlotSaveDetails('Golgi_reduced_v', project.simPlotInfo.getSimPlot('Golgi_reduced_v'), ArrayList([0]), ArrayList([0]), False, False)
    project.generatedPlotSaves.addPlotSaveDetails('Golgi_detailed_v', project.simPlotInfo.getSimPlot('Golgi_detailed_v'), ArrayList([0]), ArrayList([0]), False, False)
    # set up voltage recording at current injection site
    sim_plot_dend_reduced = project.simPlotInfo.getSimPlot('Golgi_reduced_v_d')
    sim_plot_dend_reduced.setSegmentId(str(seg_id_reduced))
    project.generatedPlotSaves.addPlotSaveDetails('Golgi_reduced_v_d', sim_plot_dend_reduced, ArrayList([0]), ArrayList([seg_id_reduced]), False, False)
    sim_plot_dend_detailed = project.simPlotInfo.getSimPlot('Golgi_detailed_v_d')
    sim_plot_dend_detailed.setSegmentId(str(seg_id_detailed))
    project.generatedPlotSaves.addPlotSaveDetails('Golgi_detailed_v_d', sim_plot_dend_detailed, ArrayList([0]), ArrayList([seg_id_detailed]), False, False)

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
