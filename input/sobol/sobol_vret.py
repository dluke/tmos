
from SALib.sample import saltelli

import sobol
import command
import numpy as np

# local configuration
import base_config

"""
fix missing
"""

reduced_names = [
	'k_ext_off', 
	'dwell_time',
	'pilivar', 
	'kb_sh',
	'k_spawn', 
	'k_ext_on',
	'k_resample',
	'k_ret_off',
	'Lp',
	'ks',
	'f_stall'
]
 
d_bound = 0.004 # nm

problem = {
	'num_vars': 11,
	'names': reduced_names,
	'bounds': [
		[0.2,1.0],
		[0.5,3.0],
		[1.0,15.0],
		[0.01/d_bound, 1.0/d_bound],
		[0.1,8.0],
		[0.20,0.50],
		[0.5, 10.0],
		[0.05, 0.20],
		[1.0, 10.0],
		[1000, 20000],
		[20, 200]
	]
}

	
def new_sampling(problem, N):
	values = saltelli.sample(problem, N, calc_second_order=False)
	sobol.write_sampled(values)


def run(settings):
	problem = sobol.read_problem()
	sampled_values = sobol.read_sampled()

	# initialize tmos wrapper
	pool = sobol.SimulationPool()
	pool.cleanup = True # delete data as we go.
	pool.problem = problem
	pool.sampled_params = sampled_values
	pool.step_target = 1000
	pool.args = base_config.args
	pool.maxcores = settings.get('maxcores', 8)
	if 'simtime' in settings:
		pool.args.system.simtime = settings['simtime']

	jobs = pool.make_jobs()
	pool.parallel_run(jobs)

if __name__ == '__main__':

	sobol.write_problem(problem)
	problem = sobol.read_problem()
	N = 2048
	new_sampling(problem, N)

	settings = {
			'simtime': 2000,
			'maxcores': 31
			}
	print("Setup complete. Attempt run...")
	run(settings)

	print("Finished")

