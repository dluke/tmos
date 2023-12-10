
import os, sys
import twutils
import numpy as np



# generator for sequential names
# NOT USED

class SeqNamer(object):
    def __init__(self, seq, pdir='plots', name='defaulttmp', extra='', pad=4):
        self.seq = iter(seq)
        base = os.path.join(pdir, name, extra)
        assert base[-1] == '/'
        if not os.path.exists(base):
            os.mkdir(base)
        self.form = base+'_{:0'+str(pad)+'d}.png'

    def next(self, step=None):
        out = self.form.format(next(self.seq)) 
        return out



