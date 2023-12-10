#!/usr/bin/env python3

# restart run from a config.txt file
import os, sys

import parameters

def parse_random_seed():
    with open('random_seeds.txt', 'r') as f:
        seeds = []
        for line in f:
            idx, seed = map(int, line.strip().split())
            seeds.append(seed)
    return seeds

def reargs():
    args = parameters.thisread().apply()
    if not os.path.exists('config.txt'):
        sys.exit("Attempted to rerun simulation but failed to find config.txt file")
    with open('save_config.txt', 'w') as f:
        f.write(str(args))

    args.set_comment("rerun")
    #args.system.vtkwrite = True
    args.system.repeats = 1
    args.system.randominit = False
    #args.system.simtime = 200
    #args.Cell3d.repulsive_only = False
    # args.system.debug = True
    return args

def singlerun():
    args = reargs()
    import runtfp_global

    # go 
    runtfp_global.single_run(args)

def pyrun(par, kw={}):
    args = reargs()
    for k, v in kw.items():
        print('modifying parameter {}={}'.format(k,v))
        args.pset(k, v)
    #
    if par.use_seed:
        seed = parse_random_seed()[0]
        print('reusing seed ', seed)
        args.system.seed = seed
    import ctfp3d
    ctfp3d.pyrun(args, None)


if __name__=='__main__':

    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('--use-seed', action='store_true', default=True)
    parse.add_argument('-t', '--simtime', type=float)
    parse.add_argument('-track', action='store_true', default=None)
    mainloop = parse.add_mutually_exclusive_group()
    mainloop.add_argument('--pyrun', action='store_true')
    mainloop.add_argument('--single', action='store_true')
    mainloop.add_argument('--run-with-vtk', action='store_true')
    par = parse.parse_args()
    print(par)

    ## 
    kw = {}
    if par.simtime:
        kw['simtime']= par.simtime
    if par.track:
        kw['track'] = par.track
    ## 
        
    if par.pyrun:
        pyrun(par, kw)
    elif par.single:
        singlerun()
    if par.run_with_vtk:
        kw['vtkwrite'] = True
        pyrun(par, kw)

