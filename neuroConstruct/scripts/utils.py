import os
import time
import subprocess
from collections import deque

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
