
from collections import OrderedDict
from glob import glob
import sys, os

import scipy.io
import numpy as np
from tqdm import tqdm


import plotutils as utils
import twutils
import readtrack


# define constants and names for operating with matfile
# we resuse the names for the analysis module
import matdef

def getn(fname):
    prefix = 'allData_sectioin'
    ext = '.mat'
    fname = os.path.basename(fname)
    return int(fname[len(prefix):-len(ext)])

"""
 Define an interface to Fanjin data.
 tracks have different lengths but always have the same fixed time resolution
 Different length tracks means we can't use numpy array blocks
 Unlike twanalyse.Track class which uses an ordered dictionary of numpy arrays
 this class manages a set of sets of tracks.
"""

# convience constructor
def this_fjtrack():
    tracks = readtrack.trackset()
    Ftracking = ExpTracks.usesimulated(tracks)
    return Ftracking

# keep track of which file the data is in aswell as which day it was produced
# should check consistency of all calculation across files
class ExpTracks(object):

    SOURCE_FJ = "FanJin"
    SOURCE_SIM = "simulated"
    SOURCE_YR = "Yow-Ren"

    npyform = 'exptrack_{:05d}.npy'

    def __init__(self):
        """
        Use numpy array blocks to store tracking data for fast slicing.
        Keep track of a dictionary which maps column master names onto array indices
        """
        self.colmap = {}
        self.tracks = []
        self.ntracks = 0
        self.timestep = 0.1 # important default for simulated data
        self.source = None
        # for FanJin experimental data
        self.datas = None
        self.datadirs = None

    ##########################################################################
    # using numpy blocks for faster loading

    def writenpy(self, npydir):
        """Write all the numpy track objects to file for faster retrieval"""
        form = os.path.join(npydir, ExpTracks.npyform)
        for i, track in enumerate(self.tracks):
            out = form.format(i)
            print("Saving to ", out)
            np.save(out, track)

    def get_total_duration(self):
        # the total duration of all the tracks end to end 
        return np.sum(track.shape[0] for track in self.tracks) * self.timestep

    @classmethod
    def load_one(cls, npydir, idx):
        file = os.path.join(npydir, ExpTracks.npyform).format(idx)
        track = np.load(file)
        new = cls()
        new.tracks = [track]
        new.ntracks = len(new.tracks)
        new.source = ExpTracks.SOURCE_FJ
        new.colmap = matdef.fjcolmap.copy()
        new.nidx = [idx]
        return new

    @classmethod
    def load_all(cls, npydir, nidx):
        tracks = []
        for idx in tqdm(nidx):
            file = os.path.join(npydir, ExpTracks.npyform).format(idx)
            track = np.load(file)
            tracks.append(track)
        new = cls()
        new.tracks = tracks
        new.ntracks = len(new.tracks)
        new.source = ExpTracks.SOURCE_FJ
        new.colmap = matdef.fjcolmap.copy()
        # TODO establish this member in init()
        new.nidx = nidx
        return new


    ##########################################################################

    # convert to simulated data track format
    def get_tracklike(self):
        # make sure we have ax_x, ax_y, trail_x, trail_y and x, y
        self.compute_orientation()
        self.compute_ax()
        self.add_zcolumn() 
        required = ['ax_x', 'ax_y', 'x', 'y', 'trail_x', 'trail_y']
        assert( np.all([name in self.colmap for name in required]) )
        # create structured array
        dt= matdef.make_dtype(list(self.colmap.keys()))
        def make_structured_array(track):
            struct = np.empty(track.shape[0], dtype=dt)
            if False:
                print('create structured array with length {} using track with shape {}'.format(len(struct), track.shape))
            for i, head in enumerate(dt.names):
                struct[head] = np.array(track[:,self.colmap[head]], dtype=dt[head])
            return struct
        #
        print("Converting tracks to simulation format")
        trs = [readtrack.TrackLike(make_structured_array(track), ExpTracks.SOURCE_FJ) for track in tqdm(self.tracks)]
        return trs

    def track_length(self):
        """All track lengths"""
        lengths = []
        for track in self.tracks:
            lengths.append(track.shape[0])
        return lengths

    def count_tracks(self):
        return len(self.tracks)

    def filter_by_method(self, toaccept):
        """toaccept is a filtering method which takes a track and returns True/False"""
        self.tracks[:] = [track for track in self.tracks if toaccept(track)]
        lost = self.ntracks - len(self.tracks)
        self.ntracks = len(self.tracks)
        print('Use filter ', toaccept.__name__)
        print('->filtered out {} tracks. Tracks remaining {}.'.format(lost, self.ntracks))

    def compute_cxy(self):
        """compute center_x, center_y from trailxy and (head) xy
        add the new columns to the track data arrays and update self.colmap
        """
        cm = self.colmap
        for i, track in enumerate(self.tracks):
            cx = (track[:,cm['trail_x']] + track[:,cm['x']])/2.
            cy = (track[:,cm['trail_y']] + track[:,cm['y']])/2.
            cx = cx[:,np.newaxis]
            cy = cy[:,np.newaxis]
            track = np.concatenate([track, cx, cy], axis = 1)
            self.tracks[i] = track
        ncols = len(self.colmap)
        self.colmap['center_x'] = ncols
        self.colmap['center_y'] = ncols + 1

    def add_empty(self, colname):
        # cheap function for adding a dummy column
        for i, track in enumerate(self.tracks):
            ntrack = np.concatenate([track, np.empty((track.shape[0],1))], axis=1)
            self.tracks[i] = ntrack
        self.colmap[colname] = len(self.colmap)

    def compute_ax(self):
        """compute the unit length axis of the bacterium body i.e. ax_x, ax_y columns
        Using orientation column
        """
        if 'orientation' not in self.colmap:
            self.add_empty('orientation')
            self.override_orientation_using_lead_trail()
        cm = self.colmap
        for i, track in enumerate(self.tracks):
            theta = self.get_col(track, 'orientation')
            ax_x = np.cos(theta)[:,np.newaxis]
            ax_y = np.sin(theta)[:,np.newaxis]
            ntrack = np.concatenate([track, ax_x, ax_y], axis = 1)
            self.tracks[i] = ntrack
        ncols = len(self.colmap)
        self.colmap['ax_x'] = ncols
        self.colmap['ax_y'] = ncols + 1

    def compute_aspect(self):
        self.add_empty('aspect')
        for i, track in enumerate(self.tracks):
            length, width = self.get_col(track, 'length'), self.get_col(track, 'width')
            aspect = length/width
            self.tracks[i][:,-1] = aspect

    def compute_orientation(self):
        self.add_empty('orientation')
        self.override_orientation_using_lead_trail()

    def add_zcolumn(self):
        for i, track in enumerate(self.tracks):
            zero = np.zeros((track.shape[0],1))
            ntrack = np.concatenate([track, zero, np.copy(zero)], axis=1)
            self.tracks[i] = ntrack
        ncols = len(self.colmap)
        self.colmap['z'] = ncols
        self.colmap['trail_z'] = ncols + 1

    def get_cols(self, name):
        return [track[:,self.colmap[name]] for track in self.tracks]

    def get_col(self, track, name):
        return track[:,self.colmap[name]]

    def get_step(self, timestep):
        fstep = float(timestep)/self.timestep
        if not fstep.is_integer():
            raise RuntimeError("Failed to use timestep {} for data with base interval {}"
                    .format(timestep, self.timestep))
        return int(round(fstep))

    def get_columns(self, cols, timestep=None, tracknum=None):
        """
        """
        tracks =self.tracks if tracknum == None else [self.tracks[i] for i in tracknum]
        """Option to choose the specific track/s to plot by order of the list"""
        timestep = self.timestep if timestep == None else timestep
        """timestep should be an exact multiple of self.timestep"""

        step = self.get_step(timestep)
        coltracks = []
        colidx = [self.colmap[col] for col in cols]
        if len(colidx) == 1:
            [colidx] = colidx # ?
        for track in tracks:
            block = track[::step,colidx]
            # maybe unnecessary
            block = np.squeeze(block)
            coltracks.append(block.T)
        return coltracks
    

    def get_whole_col(self, col):
        return np.concatenate(self.get_columns([col]))

    def get_velocity(self, tracknum=[0]):
        """return an array with x,y velocity values
        this will be used with
        https://www.nature.com/articles/ncomms8516#Sec19 code for autoregression model
        """
        # head velocity
        headxy = self.get_columns(['x', 'y'], tracknum=tracknum)
        xy = headxy[0]
        return (xy[:,1:] - xy[:,:-1]).T/self.timestep


    # accept a list of tracks, each a dictionary of columns
    @classmethod
    def usesimulated(cls, trs):
        # join the columns together and make them readable with a new name
        tvals = ['time', 'x','y', 'trail_x', 'trail_y']
        def get_stack(tr):
            columns = np.column_stack([tr.track[t] for t in tvals])
            # convert from time to simple timeid
            # I can't remeber if timeids start at zero or one in the matlab data
            size = columns.shape[0]
            columns[:,0] = np.arange(1, size+1, 1)
            return columns

        new = cls()
        new.colmap = dict(list(zip(tvals, list(range(len(tvals))))))
        new.tracks = [get_stack(tr) for tr in trs]
        new.ntracks = new.count_tracks()
        print('Create FJ track object from simulated data')
        print('the shape of the first track is', new.tracks[0].shape)
        new.source = ExpTracks.SOURCE_SIM
        new.verify()
        return new

    @classmethod
    def usefj(cls, datadirs, first=True, tracktype=matdef.DENOISED):
        """Class constructor for FanJin Data
        2. cls.loadfirst() / cls.loadall()
        3. cls.compile_tracks()
        """
        new = cls()
        new.colmap = matdef.fjcolmap
        new.datadirs = datadirs
        new.files, new.keymap = new.find_files()
        new.datas = []
        new.source = ExpTracks.SOURCE_FJ
        if first:
            new.loadfirst()
        else:
            new.loadall()
        new.compile_tracks(tracktype=tracktype)
        new.verify()
        #new.filter_smooth_orientation()
        new.override_orientation_using_lead_trail()
        return new
    
    @staticmethod
    def _yr_get_distances(colmap):
        # get column names which correspond to distances (for converting to \mu m)
        distanceidx = []
        for dname in matdef.yrdistance_names:
            if dname in colmap:
                distanceidx.append(colmap[dname])
        return np.array(distanceidx)

    @classmethod
    def useyr(cls, yrfile, colmap=matdef.yrcolmap, yrpix_to_mu=0.042, timestep=30.,
            convert=True):
        """convert=False is just for checking the columns are read correctly"""
        new = cls()
        new.timestep = timestep
        new.colmap = colmap
        yrdata = np.loadtxt(yrfile)
        if convert:
            # converted from frameid to an actual time
            yrdata[:,colmap['time']] *= timestep
            # converted pixels to distance
            yrdata[:,cls._yr_get_distances(colmap)] *= yrpix_to_mu
            # convert to radians
            yrdata[:,colmap['orientation']] *= np.pi/180.

        centries = yrdata[:,-1]
        cellids = np.unique(centries)
        new.tracks = []
        for cid in cellids:
            # slow method
            # To speed up. Efficiently find first instance of each cell id in file
            cargs = np.argwhere(centries == cid)
            new.tracks.append(np.squeeze(yrdata[cargs,:]))

        new.source = ExpTracks.SOURCE_YR
        new.ntracks = len(new.tracks)
        new.timestep = timestep
        new.verify()
        return new

    def persistent_orientation(self):
        """Determine a constistent polar orientation for YR data"""
        flip_threshold = np.pi - np.pi/20
        for trackdata in self.tracks:
            ors = trackdata[:,self.colmap['orientation']]
            postdiff = np.insert(ors[1:] - ors[:-1], 0, 0)
            """insert 0 at the start so that flip positions are on the frame following the flip"""
            postflip = np.argwhere(np.abs(postdiff) > flip_threshold)
            if postflip.size == 0:
                continue
            elif postflip.size % 2 == 1:
                postflip = np.append(postflip, None) # implicit dtype = object
            for evenpairs in postflip.reshape(-1,2)[::2]:
                fr, to = evenpairs
                correction = np.pi if postdiff[fr] <= 0 else -np.pi
                """This correction orientation to (-pi, pi] range"""
                ors[fr:to] + correction
            trackdata[:,self.colmap['orientation']] = ors

    def choose_head_orient(self):
        """
        process the orientation column again so that taking cos and sin of theta 
        gives a vector oriented towards the head.
        """

    def verify(self):
        """Verify that the necessary fields have been set for construction of the object"""
        required = [self.colmap, self.tracks, self.ntracks, self.source]
        safe = all([bool(field) for field in required])
        if not safe:
            raise RuntimeError("Failed to verify construction of ExpTracks object")

    #################################################################
    # FJ data manipulation 

    def find_files(self):
        allfiles = []
        filekeys= []
        for datadir in self.datadirs:
            path = os.path.join(datadir, '*.mat')
            sortfiles= sorted(glob(path), key=getn)
            allfiles.extend(sortfiles)
            for f in sortfiles:
                filekeys.append((datadir, getn(f)))

        keymap= OrderedDict(list(zip(filekeys, list(range(len(allfiles))))))
        return allfiles, keymap

    def load(self, fkey):
        matfile = self.files[self.keymap[fkey]]
        print('loading data from file {}'.format(matfile))
        thisdata = scipy.io.loadmat(matfile)
        self.datas.append(thisdata)

    def loadfirst(self):
        self.load(list(self.keymap.keys())[0])
    # now operate on what we chose to load
    def loadall(self):
        for fkey in self.keymap:
            self.load(fkey)

    @staticmethod
    def data_tracks(data):
        topkey = 'allData_section'
        return data[topkey][0][0][0][0]

    def override_orientation_using_lead_trail(self):
        lts = self.get_columns(matdef.PLEADTRAIL)
        for i, lt in enumerate(lts):
            time, x, y, tx, ty = lt
            xd = x - tx
            yd = y - ty
            orientation = np.arctan2(yd,xd)
            self.tracks[i][:,self.colmap['orientation']] = orientation


    def compile_tracks(self, tracktype=matdef.DENOISED):
        # read in all the data to self.tracks which is list of arrays
        for data in self.datas:
            data_tracks = self.data_tracks(data)
            print('.mat file has {} tracks'.format(len(data_tracks)))
            for dtrack in data_tracks:
                # select the denoised data or original data
                track = dtrack[0][0][tracktype]
                # convert from a timeid to an actual timestep
                track[:,self.colmap['time']] *= matdef.TIMESTEP
                track[:,self.colmap['orientation']] *= np.pi/180
                self.tracks.append(track)
        self.ntracks = len(self.tracks)


    #################################################################################
    # filters
    def filter_smooth_orientation(self, ismooth=np.pi/40.):
        """
        Atleast one of the tracks appears as a stationary circle when plotted
        and has badly defined orientation.
        Define a criteria to exlude such tracks"""
        def accept_smooth(track):
            ors = track[:,self.colmap['orientation']]
            dors = np.mean(np.abs(ors[1:] - ors[:-1]))
            return dors < ismooth
        self.filter_by_method(accept_smooth)

    def filter_short(self, mtime):
        mstep = self.get_step(mtime)
        def toaccept(track):
            return True if track.shape[0] >= mstep else False
        self.filter_by_method(toaccept)


