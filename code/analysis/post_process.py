#!/usr/bin/env python3

import wrinterval as wr
from readtrack import find_track_files

if __name__=='__main__':
    ptr_files = find_track_files(single='pilitracking.dat', form='pili_*.dat', ddir='data/')
    wr.post_process(ptr_files, wr._bound_condition, wr.nameform)
