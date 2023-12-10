
import sobol
import command
import numpy as np

# local configuration
import base_config

"""
 lets assume k_ret_off is large enough that pili always complete a full retraction
 and that retraction is determined either by surface interaction or by k_ext_off
 k_ext_on is also irrelevant in the one retraction per pilus limit
 that leaves k_ext_off as the only important parameter and it is relevant only in determining
 attachments of long pili (should exactly determine the average length of non interacting pili)
"""

reduced_names = [
    'dwell_time',
    'pilivar', 
    'kb_sh',
    'k_spawn', 
]

d_bound = 0.004 # \mu m
 
reduced_problem = {
    'num_vars': 4,
    'names': reduced_names,
    'bounds': [
        [0.05,3.0],
        [1.0,15.0],
        [0.01/d_bound, 1.0/d_bound],
        [0.1,8.0]
    ]
}
    
def new_sampling(problem, N):
    values = sobol.mcsampling(problem, N)
    sobol.write_sampled(values)

def run(settings):
    problem = sobol.read_problem()
    sampled_values = sobol.read_sampled()

    # initialize tmos wrapper
    pool = sobol.SimulationPool()
    pool.cleanup = True
    pool.problem = problem
    pool.step_target = 1000
    pool.sampled_params = sampled_values
    pool.args = base_config.args
    pool.maxcores = settings.get('maxcores', 8)
    if 'simtime' in settings:
        pool.args.system.simtime = settings['simtime']

    jobs = pool.make_jobs()
    pool.parallel_run(jobs)


if __name__ == '__main__':

    sobol.write_problem(reduced_problem)
    problem = sobol.read_problem()
    N = 10000
    new_sampling(problem, N)

    settings = {
            'simtime': 2000,
            'maxcores': 64
            }
    print("Setup complete. Attempt run...")
    run(settings)

