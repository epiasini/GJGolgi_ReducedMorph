import os
import time
import subprocess
from collections import deque

def pull_remotes(sim_refs, general_sim_dir='../simulations/'):
    for sim_ref in sim_refs:
	sim_path = general_sim_dir + sim_ref
	pullsimfile_path = sim_path + '/pullsim.sh'
	print('Pulling from ' + sim_ref)
	subprocess.call([pullsimfile_path])

def wait_and_pull_remote(sim_refs, sleep_time=5, general_sim_dir='../simulations/'):
    sim_refs = deque(sim_refs)
    while sim_refs:
        # keep on going through the queue of running sims, and try to pull
	# the results from the cluster. sim_refs can a list as well as a deque.
	sim_ref = sim_refs.popleft()
	sim_path = general_sim_dir + sim_ref
	pullsimfile_path = sim_path + '/pullsim.sh'
	timefile_path = sim_path + '/time.dat'
	print('Pulling from ' + sim_ref)
	subprocess.call([pullsimfile_path])
	if not os.path.exists(timefile_path):
	    sim_refs.append(sim_ref)
	time.sleep(sleep_time)

def ir_sim_ref(timestamp, gj_conn_type, amplitude, trial):
    return 'ir' + '_' + timestamp + '_' + gj_conn_type + '_' + str(int(round(amplitude))) + '_t' + str(trial)

def ir_single_cell_sim_ref(timestamp, amplitude):
    return 'ir_sc' + '_' + timestamp + '_' + str(int(round(amplitude)))

def cs_sim_ref(timestamp, gj_conn_type, cell, trial):
    return 'cs' + '_' + timestamp + '_' + gj_conn_type + '_c' + str(cell) + '_t' + str(trial)

def cs_cell_positions_file(timestamp, gj_conn_type, trial):
    return '../scriptedSimulations/cs' + '_' + timestamp + '_' + gj_conn_type + '_t' + str(trial) + '_positions.csv'  

def cs_edge_list_file(timestamp, gj_conn_type, trial):
    return '../scriptedSimulations/cs' + '_' + timestamp + '_' + gj_conn_type + '_t' + str(trial) + '_edges.csv'  

def coupling_coefficient(r01, rl0, rl1, dv, I):
    return (rl1/(rl1+r01)) - dv*r01*(rl0+rl1+r01) / ((rl1+r01) * (I*rl0*(rl1+r01) + dv*rl0))

def variable_heterogeneity(timestamp, deg_mean, trial):
    return 'vh' + '_' + timestamp + '_m' + str(deg_mean) + '_t' + str(trial)
