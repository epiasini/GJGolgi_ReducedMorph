import os
import time
import random
import networkx as nx
from math import fabs

from java.awt import Color
from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol import neuroconstruct as nc

import utils

sim_duration = 2000

timestamp = str(time.time())
pm = nc.project.ProjectManager(None,None)
project_path = '../GJGolgi_ReducedMorph.ncx'
project_file = File(project_path)
project = pm.loadProject(project_file)

remote_sim_refs = []

for year in ['2010', '2012']:

    sim_config_name = 'network_desynchronisation_reduced_only_' + year + 'gap'
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
    sim_config.setSimDuration(sim_duration)
    project.neuronSettings.setNoConsole()

    # generate
    nC_seed = random.getrandbits(32)
    pm.doGenerate(sim_config_name, nC_seed)
    while pm.isGenerating():
        time.sleep(0.02)
    print('network generated')

    sim_ref = 'rd_' + timestamp + '_y' + year
    project.simulationParameters.setReference(sim_ref)

    # generate and compile neuron files
    print "Generating NEURON scripts..."
    project.neuronFileManager.setSuggestedRemoteRunTime(40)
    simulator_seed = random.getrandbits(32)
    project.neuronFileManager.generateTheNeuronFiles(sim_config,
                                                     None,
                                                     nc.neuron.NeuronFileManager.RUN_HOC,
                                                     simulator_seed)

    compile_process = nc.nmodleditor.processes.ProcessManager(project.neuronFileManager.getMainHocFile())
    compile_success = compile_process.compileFileWithNeuron(0,0)


    # simulate
    if compile_success:
        print "Submitting simulation reference " + sim_ref
        pm.doRunNeuron(sim_config)
        time.sleep(2.5) # Wait for sim to be kicked off
        if sim_config.getMpiConf().isRemotelyExecuted():
            remote_sim_refs.append(sim_ref)
        else:
            # if running locally, never have more than one sim running
            # at the same time
            print('Simulating on the local machine.')
            timefile_path = '../simulations/' + sim_ref + '/time.dat'
            while not os.path.exists(timefile_path):
                time.sleep(0.5)

if remote_sim_refs:
    utils.wait_and_pull_remote(remote_sim_refs, sleep_time=0.5)  

print('timestamp ' + timestamp)
System.exit(0)

