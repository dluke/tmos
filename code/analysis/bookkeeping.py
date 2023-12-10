
import os
from glob import glob
import json

import stats

# we want to know

# simulation time, execution time (clock time), data usage (breakdown according to each output type)


# example timing file
"""
cell 0000 simtime   2000.0084 clock_time 31.41074776649475
cell 0001 simtime   2000.5446 clock_time 37.75506615638733
cell 0002 simtime   2000.0032 clock_time 34.85314726829529
cell 0003 simtime   2012.6875 clock_time 33.194913387298584
cell 0004 simtime   2000.0117 clock_time 32.31716966629028
cell 0005 simtime   2000.5105 clock_time 27.00812005996704
cell 0006 simtime   2005.3936 clock_time 31.604656219482422
cell 0007 simtime   2000.1407 clock_time 26.362374782562256
cell 0008 simtime   2022.7877 clock_time 35.02218008041382
cell 0009 simtime   2000.3109 clock_time 34.23537302017212
"""

# assume standard names
trdataform = 'bacterium_*.dat'
ptrdataform = 'pili_*.dat'
etrdataform = 'event_*.dat'
datatypes = ['tracking', 'pili', 'event']

# use stats
@stats.keep
def summary():
    cellid, simtimes, clocktimes = read_timing()
    sd = {}
    sd['total_simtime'] = sum(simtimes)
    sd['total_clocktime'] = sum(clocktimes)
    trdatasize, ptrdatasize, etrdatasize = get_data_usage('data/')
    # leave this in bytes
    sd['trdatasize'] = trdatasize
    sd['ptrdatasize'] = ptrdatasize
    sd['etrdatasize'] = etrdatasize
    sd['total_data_size'] = sum([trdatasize, ptrdatasize, etrdatasize])
    return sd

##################################################################################
# methods for timing and data usage

timingfile = 'timing.txt'

def read_timing():
    cell = []
    simtimes = []
    clocktimes = []
    type_f = [int, float, float]
    with open(timingfile, 'r') as f:
        for line in f:
            ls = line.strip().split()
            assert(ls[0] == 'cell')
            cellid, _simtime, _clocktime = [type_f[i](x) for i, x in enumerate(ls[1::2])]
            cell.append(cellid)
            simtimes.append(_simtime)
            clocktimes.append(_clocktime)

    return cellid, simtimes, clocktimes

def clocktime():
    cellid, simtimes, clocktimes = read_timing()
    return sum(clocktimes)

def simtime():
    cellid, simtimes, clocktimes = read_timing()
    return sum(simtimes)

# no recursive get size
def _data_size(ddir):
    return sum(os.path.getsize(f) for f in os.listdir(ddir) if os.path.isfile(f))

def get_data_usage(ddir):
    trdatapaths = glob(os.path.join(ddir, trdataform))
    ptrdatapaths = glob(os.path.join(ddir, ptrdataform))
    etrdatapaths = glob(os.path.join(ddir, etrdataform))
    trdatasize = sum(os.path.getsize(f) for f in trdatapaths)
    ptrdatasize = sum(os.path.getsize(f) for f in ptrdatapaths)
    etrdatasize = sum(os.path.getsize(f) for f in etrdatapaths)
    return trdatasize, ptrdatasize, etrdatasize

def _to_gb(v):
    return float(v)/10**9



##################################################################################

def test():
    val = clocktime()
    print('total clocktime {} seconds'.format(val))

def test_get_data_usage():
    gb_form = '{:5.3f}GB'
    def pretty(val):
        return gb_form.format(_to_gb(val))
    output = [pretty(v) for v in get_data_usage('data/')]
    outdct = {k:v for k, v in zip(datatypes, output)}
    print(json.dumps(outdct, indent=4))


if __name__=='__main__':

    # test()
    # test_get_data_usage()

    summary()