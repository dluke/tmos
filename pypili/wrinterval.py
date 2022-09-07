# pili.py contains the trackxy object but this is only good for tracking one cell
# We want to just output the xy position and orientation states of an array of cells efficiently

import numpy as np
import os, sys
from glob import glob

from collections import OrderedDict, namedtuple

import txtdata

Simstate = namedtuple('Simstate', ['kmcstep', 'cell', 'istep', 'force'])

# global data directory name
datdir = 'data/'

mintrack_2d = ['time', 'x', 'y', 'ax_x', 'ax_y']
track_2d = [
        'time', 'process', 'x', 'y', 'trailx', 'traily', 'ax_x', 'ax_y',
        'nbound', 'pl_avg', 'nsteps', '|force|', 'en_stretch'
        ]

# combine the binding object here by generalising hform 
class Trackxy(object):
    def __init__(self, fname='trackxy.dat', header=None, separator='  '):
        assert header, 'Trackxy Initialisation Failed. No header'

        self.fo = open(fname, 'w')
        self.header = header
        self.separator = separator
        # define lists of format strings
        self.linelist = [txtdata.hform[head] for head in self.header]
            # write the header
        headform = [txtdata.hwidth[head] for head in self.header]
        self.fo.write(self.separator.join(headform).format(*header) + '\n')

    def what_file():
        return self.fo.name

    @classmethod
    def setup_mintrack(cls, fname):
        return cls(fname, mintrack_2d)

    def compile_line(self,lvals):
        line = self.separator.join(self.linelist) + '\n'
        linestr = line.format(*lvals)
        return linestr

    def write_line(self,lvals):
        linestr = self.compile_line(lvals)
        self.fo.write(linestr)
        return linestr

    # the minimal quantities
    def mintrack(self, time, cell):
        ax = cell.axis
        #hpt = cell.real_centre()
        hpt = cell.hexpt
        lvals = [time] + [hpt.x, hpt.y] + [ax.x, ax.y]
        self.write_line(lvals)

    # to generalise tracking, bind a whole set of getter methods which operate on simstate
    def track(self, Simstate):
        kmcstep, cell, istep, force = Simstate
        tme, process, _ , _ , _ = kmcstep
        x,y = cell.hexpt.x, cell.hexpt.y
        trail = cell.get_trail()

        nbound = cell.nbound()
        navg = cell.pl_avg()
        en_s = cell.energy()
        lvals = [tme, process, x, y, 
                trail.x, trail.y, 
                cell.axis.x, cell.axis.y, 
                nbound, navg, istep, force, en_s ] 

        return self.write_line(lvals)

    def close(self):
        self.fo.close()


attachment_event = [
    'time', 'pidx', 'process', 'statetime', 
    'pleq', 'plength', 'inter_excess', 'nseg', 'lasta', 
    'anchor_x', 'anchor_y', 'anchor_z',
    'paxis_x', 'paxis_y', 'paxis_z',
    'attach_x', 'attach_y', 'attach_z']

event_3d = [
    'time', 'pidx', 'process', 'trigger', 'inter_excess', '|force|',
    'pleq', 'plength', 'nseg', 'lasta', 'isbound', 'istaut', 'ext_motor', 'ret_motor',
    'last_shortening', 
    'anchor_x', 'anchor_y', 'anchor_z',
    'paxis_x', 'paxis_y', 'paxis_z',
    'attach_x', 'attach_y', 'attach_z']

# Subclass and override the track method for a new type of output
class Track3d(Trackxy):
    def __init__(self, fname='trackxy.dat', header=None):
        super(Track3d, self).__init__(fname, header=header)

    track_3d = [
        'time', 'process', 'pidx', 'x', 'y', 'z', 'trail_x', 'trail_y', 'trail_z',
        'ax_x', 'ax_y', 'ax_z', 'rax_x', 'rax_y', 'rax_z',
        'npili', 'nbound', 'ntaut', 'ncontacts', 'l_total', 'fluor_ntfp', 'iscat_ntfp'
        ]

        
    @classmethod
    def setup_track_3d(cls, fname):
        return cls(fname, cls.track_3d)
        
    # override the track method for my specific 3d output
    def track(self, Simstate):
        kmcstep, cell, istep, force = Simstate
        tme, process, _ , pidx = kmcstep
        r = cell.get_headpt()
        t = cell.get_trail()
        ax = cell.get_axis()
        rax = cell.body.frame.e1
        nbound = cell.nbound()
        # pl_avg = cell.pl_avg()
        l_total = cell.l_total()
        nc = cell.get_num_contacts()
        lvals = [tme, process, pidx,
                r.x, r.y, r.z,
                t.x, t.y, t.z, 
                ax.x, ax.y, ax.z,
                rax.x, rax.y, rax.z,
                len(cell.pili), cell.nbound(), cell.ntaut(), 
                nc, l_total, 
                cell.num_pili(0.3), cell.num_pili(1.0)
            ]
                # istep, force,
        return self.write_line(lvals)

    def track_attachment(self, ev):
        tme = ev.get_time()
        process = ev.process
        cell = ev.get_data()
        cell.update_anchors() # update axisEqs
        statetime = -1

        pilus = cell.get_pilus(ev.pidx)
        nseg = pilus.n if hasattr(pilus, 'n') else 0
        lasta = pilus.lasta if hasattr(pilus, 'lasta') else 0.
        # anchor is in body frame
        anchor = pilus.anchor
        # attach is in lab frame
        attach = pilus.attachment 

        vidx = cell.get_pilus_vidx(pilus)
        paxis = pilus.axisEq
        # add the surface intersection length
        inter_excess = pilus.get_last_inter_len()

        lvals = [tme, pilus.idx, process, statetime, 
                pilus.leq, pilus.pl, inter_excess, nseg, lasta, 
                anchor.x, anchor.y, anchor.z,
                paxis.x, paxis.y, paxis.z,
                attach.x, attach.y, attach.z]
        return self.write_line(lvals) 


        
    @classmethod
    def setup_event_3d(cls, fname):
        return cls(fname, event_3d)

    # Events: attach, release, istaut, isslack, dissolve, spawn
    def track_event(self, ev):
        tme = ev.get_time()
        process = ev.process
        pidx = ev.pidx
        cell = ev.get_data()
        cell.update_anchors()

        pilus = cell.get_pilus(pidx)

        anchor = pilus.anchor
        attach = pilus.attachment
        # vidx = cell.get_pilus_vidx(pilus)
        paxis = pilus.axisEq
        trigger = ev.trigger if len(ev.trigger) else 'None'
        fmag = pilus.force_l()
        nseg = pilus.n if hasattr(pilus, 'n') else 0
        lasta = pilus.lasta if hasattr(pilus, 'lasta') else 0.

        has_iexcess = trigger == 'attach' or process == 'release'
        inter_excess = pilus.get_last_inter_len() if has_iexcess else -1
        
        lvals = [tme, pidx, process, trigger, inter_excess, fmag,
                pilus.leq, pilus.pl, nseg, lasta, 
                pilus.isbound, pilus.istaut(), pilus.ext_motor, pilus.ret_motor,
                pilus.last_shortening,
                anchor.x, anchor.y, anchor.z,
                paxis.x, paxis.y, paxis.z,
                attach.x, attach.y, attach.z]
        return self.write_line(lvals)

    mdevent = [
        'time', 'pidx', 'process', 'nbound', 'pbrf', 'nsteps', 'rms', 'd_x', 'd_y', 'd_z', 'd_tx', 'd_ty', 'd_tz', 
         'd_ax_x', 'd_ax_y', 'd_ax_z', 'd_rax_x', 'd_rax_y', 'd_rax_z'
    ]

    @classmethod
    def setup_mdevent(cls, fname):
        return cls(fname, cls.mdevent)

    def track_mdevent(self, data):
        time, pidx, process, cell, summary = data
        # print(data)
        # print(summary.prior_pilus)
        # print(summary.accepted_pilus)
        # print()
        #
        prior = summary.prior_state
        accepted = summary.accepted_state
        nbound = cell.nbound()
        pbrf = cell.pbrf()
        dr = accepted.get_headpt() - prior.get_headpt()
        dtr = accepted.get_endpt() - prior.get_endpt()
        dcxy = accepted.frame.origin - prior.frame.origin
        d_ax = accepted.frame.e3 - prior.frame.e3
        d_rax = accepted.frame.e3 - prior.frame.e3
        #
        lvals = [time, pidx, process, nbound, pbrf,
                summary.nsteps, summary.rms,
                dr.x, dr.y, dr.z,
                dtr.x, dtr.y, dtr.z,
                d_ax.x, d_ax.y, d_ax.z,
                d_rax.x, d_rax.y, d_rax.z]
        return self.write_line(lvals)


pilitracking = [
        'pidx', 'ext_motor', 'ret_motor', 'isbound', 'istaut', 'cycles', 'statetime', '|force|',
        'pleq', 'plength', 'nseg', 'lasta', 
        'anchor_x', 'anchor_y', 'anchor_z',
        'paxis_x', 'paxis_y', 'paxis_z',
        'attach_x', 'attach_y', 'attach_z']

class PiliTrack3d(Trackxy):

    timeline_form = '%s  %s\n' % (txtdata.hform['time'], txtdata.hform['npili'])

    def __init__(self, fname='pilitracking.dat', header=pilitracking, include_unbound=False):
        super(PiliTrack3d, self).__init__(fname, header=header)
        self.fo.write('#Each block is headed by the time and number of pili\n')
        self.include_unbound = include_unbound
   
    def write_context(self, context):
        time, npili = context
        self.fo.write(PiliTrack3d.timeline_form.format(time, npili))

    def track_pili(self, tme, cell):
        cell.update_anchors() # update axisEqs
        self.write_context((tme, len(cell.pili)))
        for i, pilus in enumerate(cell.pili):
            if not self.include_unbound and not pilus.isbound:
                continue
            anchor = pilus.anchor
            attach = pilus.attachment
            paxis = pilus.axisEq
            nseg = pilus.n if hasattr(pilus, 'n') else 0
            lasta = pilus.lasta if hasattr(pilus, 'lasta') else 0.
            statetime = -1 # TODO 
            self.write_line([
                pilus.idx, pilus.ext_motor, pilus.ret_motor, 
                pilus.isbound, pilus.istaut(), pilus.cycles, statetime, pilus.force_l(), 
                pilus.leq, pilus.pl, nseg, lasta,
                anchor.x, anchor.y, anchor.z,
                paxis.x, paxis.y, paxis.z,
                attach.x, attach.y, attach.z
                ])
        self.fo.write('\n')


map_typecodes = {
    str : 'U16',
    int : 'i4',
    float : 'f8'
}
def make_dtype(trackfile):
    with open(trackfile, 'r') as ft:
        header = ft.readline().split()
    types = [map_typecodes[txtdata.htype[i]] for i in header]
    dtype = { 'names': header, 'formats' : types }
    return np.dtype(dtype)

#reader for the tracking data
def readtrack(trackfile):
    def hasheader():
        with open(trackfile, 'r') as ft:
            return ft.readline().split()[0].strip() == 'time'
    assert(hasheader())

    def read_block(ft, header):
        thistype = [txtdata.htype[i] for i in header]
        tracklist = [[] for _ in header]
        for line in ft:
            if len(line.split()) != len(header):
                print('WARNING: data file {} has malformed line of length {}/{}'.format(
                    ft.name, len(line.split()), len(header)))
                continue
            lvals = [thistype[i](lv) for i, lv in enumerate(line.split())]
            for i, lv in enumerate(lvals):
                tracklist[i].append(lv)
        return tracklist
        
    with open(trackfile, 'r') as ft:
        header = ft.readline().strip().split()
        tracklist = read_block(ft, header)

    # convert to structured array
    dt = make_dtype(trackfile)
    track = np.empty(len(tracklist[0]), dtype=dt)
    for i, head in enumerate(header):
        track[head] = np.array(tracklist[i], dtype=dt[head])
    return track


# reader for pili data (pili.dat files)
# returns a list [(Context, Block)]
# where context is (time, npili)
# and block is column data and can be accessed like block['pidx']
def readpilitrack(pilitrackfile):
    if pilitrackfile.endswith('.dat'):
        return _readpilitrack(pilitrackfile)
    elif pilitrackfile.endswith('.pkl'):
        return _readpkltrack(pilitrackfile)
    else:
        raise RuntimeError("Unrecognised filetype {}".format(pilitrackfile))

def _readpilitrack(pilitrackfile):
    def hasheader():
        with open(pilitrackfile, 'r') as ft:
            return ft.readline().split()[0].strip() == 'pidx'
    assert(hasheader())
    dt = make_dtype(pilitrackfile)

    def read_block(ft, header):
        thistype = [txtdata.htype[i] for i in header]
        tracklist = [[] for _ in header]
        for line in ft:
            if line.strip() == '':
                break
            lvals = [thistype[i](lv) for i, lv in enumerate(line.split())]
            for i, lv in enumerate(lvals):
                tracklist[i].append(lv)
        # convert to numpy arrays
        return [np.array(col) for col in tracklist]

    pilitrack = []

    with open(pilitrackfile, 'r') as ft:
        line = iter(ft.readlines())
        header = next(line).strip().split()
        comment = next(line)

        head = next(line)
        while head:
            time, npili = head.strip().split()
            context = (float(time), int(npili))
            block = read_block(line, header)
            # convert block to structure array (?)
            nameblock = np.empty(block[0].size, dtype=dt)
            for i, col in enumerate(block):
                try:
                    nameblock[header[i]] = col
                except:
                    print(col)
                    raise
            pilitrack.append( (context, nameblock) )
            # setup for the next loop
            try:
                head = next(line)
            except StopIteration:
                head = None

    return pilitrack

def _readpkltrack(pilitrackfile):
    assert(pilitrackfile.endswith('.pkl'))
    assert(os.path.dirname(pilitrackfile).endswith('post'))
    with open(pilitrackfile, 'rb') as fp:
        ptdata = pickle.load(fp)
    return ptdata

# post_processing pili tracking data
# (i) apply a condition to keep the data size in memory down 

def _bound_condition(block):
    return block[block['isbound'] == 1]

def pilitrack_condition(ptdata, condition):
    newdata = []
    for item in ptdata:
        context, block = item
        time, npili = context
        # choice here on whether to leave in empty blocks, I leave them in ...
        block = condition(block)
        # should I update the context as well?
        npili = len(block)
        context = (time, npili)
        newdata.append((context, block))
    return newdata

#  want to write back the processed version of the data so its faster to read
import pickle
post_process_dir = 'post/'
nameform = 'bound_pili_{:05d}.pkl'

def post_process(ptr_files, condition, nameform):
    for i, target in enumerate(ptr_files):
        ptdata = readpilitrack(target)
        newdata = pilitrack_condition(ptdata, condition)
        datadir = os.path.dirname(target)
        postdir = os.path.join(datadir, post_process_dir)
        if not os.path.exists(postdir):
            os.makedirs(postdir)
        out_form = os.path.join(postdir, nameform)
        with open(out_form.format(i), 'wb') as fp:
            print('dumping processed output to ', out_form.format(i))
            pickle.dump(newdata, fp)



###############################################################################

