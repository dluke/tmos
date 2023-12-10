
import sys

import ctfp3d as tfp
from pypili_on_plane import *

import pypili_on_plane as inputfile

args.system.simtime = 10
args.system.repeats = 1


if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--single', default=False, action='store_true', 
                        help='run the unmodified parameter config')
    parser.add_argument('--pyrun', default=False, action='store_true', 
                        help='run the unmodified parameter config with pyrun loop')
    parser.add_argument('--first', default=False, action='store_true', 
                        help='run only the first parameter set')
    parser.add_argument('--short', default=False, action='store_true', 
                        help='run a short simulation for every parameter config')

    command = parser.parse_args()
    if command.single:
        inputfile.single_run()
        sys.exit()
    elif command.pyrun:
        tfp.pyrun(args, None)
        sys.exit()

    if command.short:
        pool()
    else:
        pool(dry_run=True)

