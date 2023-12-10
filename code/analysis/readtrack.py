
# need a separate file to store the Track class for reading and operating
# on simulated data.

import sys, os
import numpy as np
import pandas as pd
import filesystem
norm = np.linalg.norm

from collections import OrderedDict
from glob import glob

import wrinterval as wr

verbose = False

#######################################################j
# i/o methods

def read_last_track(ddir):
	fpath = os.path.join(ddir, 'data/', 'bacterium_*.dat')
	last_path = sorted(glob(fpath))[-1]
	return  Track(last_path)

def get_track_file(idx, ddir = ''):
	#print os.path.join(ddir, 'bacterium_*.dat')
	trackfs = sorted(glob(os.path.join(ddir, 'bacterium_*.dat')))
	if len(trackfs) > 0:
		return trackfs[idx]
	return []

# get the filename list
# the cut parameter lets us choose a subset of the tracks under ./data
# return the filenames
def find_track_files(single='trackxy.dat', form='bacterium_*.dat', cut=(0,None), ddir=None):
	if ddir is None:
		ddir = wr.datdir # using './data'
	_form = os.path.join(ddir, form)
	if os.path.exists(single):
		print('found {} file'.format(single))
		# look for trackxy.dat
		return [single]
	elif os.path.exists(ddir):
		if verbose:
			print('searching for tracks with form {}'.format(_form))
		trackfs = sorted(glob(_form))
		if not trackfs:
			_form = os.path.join(ddir, "data/", form)
			trackfs = sorted(glob(_form))
		# look under data/
		fr, to = cut
		if to != None and to > len(trackfs):
			print('Warning: asked for {} tracks but found only {}'.format(n, len(trackfs)))
			n = None # slicing with None returns the whole sequence
		return trackfs[fr:to]
	else:
		# no data found
		print('no data found at ', os.path.abspath(ddir))
		return None

def trackset(ddir=None):
	# ddir is the data/ directory
	"""Return a list of Track objects. The complete tracking data"""
	def read_and_process(tf):
		# return Track(tf).extend_projected_axis_by_radius()
		return Track(tf)
	trs = [read_and_process(tf) for tf in find_track_files(form='bacterium_*.dat', ddir=ddir)]
	return trs

def eventset(ddir=None):
	"""Return a list of Track objects. The complete tracking data"""
	trs = [Track(tf) for tf in find_track_files(single='event.dat', form='event_*.dat', ddir=ddir)]
	return trs

def mdeventset(ddir=None):
	"""Return a list of Track objects. The complete tracking data"""
	# single mdevent file exists?
	trs = [Track(tf) for tf in find_track_files(single='mdevent.dat', form='mdevent_*.dat', ddir=ddir)]
	return trs

# def attset():
	# return [Track(os.path.join(d,'attachment.dat')) for d in filesystem.getdirs()]

def piliset(n=None, ddir=None):
	found = find_track_files(single='pilitracking.dat', form='pili_*.dat', ddir=ddir)
	return [wr.readpilitrack(tf) for tf in found]

# read back processed and pickled data
def bound_piliset():
	form = 'bound_pili_*.pkl'
	found = find_track_files(form=form, ddir='data/post/')
	print('found {} post_processed tracks'.format(len(found)))
	return [wr.readpilitrack(tf) for tf in found]

# class which understands how to slice the track 
class Track(object):

	verbose = False

	def _find_trackfile(self):
		if os.path.exists('trackxy.dat'):
			trackfile = 'trackxy.dat' 
		else:
			trackfile = get_track_file(0, wr.datdir)
		return trackfile

	def __init__(self, trackfile=None):
		self.trackfile = trackfile if trackfile else self._find_trackfile()
		
		if self.verbose:
			print('reading track from {}'.format(
					os.path.join(filesystem.thisdir(), self.trackfile)))

		# read structured array
		self._track = wr.readtrack(self.trackfile)
		self.source = "simulated"
		self._init()
		self.step_idx = None

	def _init(self):
		time = self._track['time']
		if time.size > 0:
			self.timebounds = time[0], time[-1]
			self.tstep = np.mean(time[1:] - time[:-1])
		# either a slice object or an numpy index array
		self.slice = slice(None)

	def short(self, N):
		return self.cut(0, N)

	def cut(self, a, b):
		new = self.copy()
		new._track  = self._track[a:b]
		return new

	def cut_time(self, a, b):
		dt = self['time'][1] - self['time'][0]
		ta, tb = int(a/dt), int(b/dt)
		return self.cut(ta, tb)

	def get_n2(self):
		return np.stack([self['x'], self['y']]).T


	

	# force self.track['x'] and self['x'] to be the same for backwards compatibility
	@property 
	def track(self):
		return self

	@property 
	def size(self):
		return self['time'].size

	@property 
	def x(self):
		return self['x']

	@property 
	def y(self):
		return self['y']

	def __len__(self):
		return self.size

	def iter(self):
		return self._track

	def get_structured_array(self):
		return self._track

	def get_dtype(self):
		return self._track.dtype

	def __getitem__(self, k):
		try:
			return self._track[k][self.slice]
		except IndexError:
			print(self.slice)
			raise
		
	def __setitem__(self, k, v):
		try:
			self._track[k][self.slice] = v
		except IndexError:
			print(self.slice)
			raise

	def get_dt(self):
		# actually should be careful with this ... approximate timestep
		return self['time'][1:] - self['time'][:-1]

	def get_duration(self):
		return self['time'][-1] - self['time'][0]
	
	def _clean_bad_dt(self):

		dt = self.get_dt()
		idx = (dt < 0.05)
		self._track = self._track[:][1:][np.array(range(self._track.size-1))[~idx]]

	def cut_start(self, initial_time=10): # in seconds
		in_idx = np.searchsorted(self._track['time'], initial_time)
		if in_idx == self._track.size-1:
			print("WARN: simulation is < initial cut time. refusing to cut anything")
		self._track = self._track[in_idx:]
		return self

	# filting data with a set time resolution
	# WARNING: resolution should be >> resolution of the data, data should have regular output
	# WARNING: dont use this method for filtering args.system.track = True (full output)
	def filter_to_resolution(self, resolution, initial=0):
		# honestly why did I implement this as a regular interval
		if not (initial < self._track.size):
			raise ValueError("attempting to cut {} steps from time series of size {}".format(initial, self._track.size))
		assert(resolution > self.tstep)
		time = self._track['time']
		nvals = int(round((time[-1] - time[initial])/resolution))
		args = np.array(np.linspace(initial, time.size-1, nvals+1, True), dtype=int)
		self.slice = args

	def true_filter_to_resolution(self, resolution, initial=10):
		time = self._track['time']
		# cutting off the start is important because time between steps can be large with no pili
		startidx = np.searchsorted(time, initial, side='right')
		_slice = [startidx]
		target = iter(np.arange(initial+resolution, time[-1], resolution))
		current = next(target)
		for i, t in enumerate(time[startidx:]):
			if t > current:
				_slice.append(startidx + i)
				try:
					current = next(target)
				except StopIteration:
					break
		self.slice = np.array(_slice)

	def clear_filter(self):
		self.slice = slice(None)


	def get_slice_length(self):
		return self._track[:][self.slice].size
	
	def weighted_average(self, col):
		time = self['time']
		colval = self[col]
		return np.sum(time*colval)/np.sum(time)

	def step_to_resolution(self, wsize):
		self.slice = slice(0,None,wsize)

	def filter_by_process(self, process):
		procs = np.argwhere(self['process'] == process).flatten()
		self.slice = procs

	def part(self, start=0, stop=None):
		tr = self.copy()
		tr._track = self._track[start:stop]
		return tr

	def flip_poles(self):
		self['x'], self['trail_x'] = self['trail_x'], self['x'].copy()
		self['y'], self['trail_y'] = self['trail_y'], self['y'].copy()
		self['z'], self['trail_z'] = self['trail_z'], self['z'].copy()
		self.step_idx, self.trail_step_idx = self.trail_step_idx, self.step_idx

	def get_dx(self):
		position =  self.get_head2d()
		return position[1:] - position[:-1]

	def get_body_projection(self):
		apole = self.get_head()[:,:2]
		bpole = self.get_trail()[:,:2]
		body = apole - bpole
		return body
	

	#################################################################################
	# assume linearised track, (step_idx is not None)

	def get_nsteps(self):
		return len(self.step_idx) - 1

	def get_step_dt(self, trail=False):
		step_idx = self.step_idx if not trail else self.trail_step_idx
		s_time = self['time'][step_idx]
		return s_time[1:] - s_time[:-1]
	
	def get_step_dx(self, trail=False):
		a, b = ('x', 'y') if not trail else ('trail_x', 'trail_y')
		step_idx = self.step_idx if not trail else self.trail_step_idx
		xy = np.stack([self[a][step_idx], self[b][step_idx]], axis=1)
		dx = xy[1:] - xy[:-1]
		return dx

	def get_step_velocity(self, trail=False):
		dt = self.get_step_dt(trail)
		dx = self.get_step_dx(trail)
		return np.nan_to_num( dx / dt[:,np.newaxis] )

	def get_step_speed(self, trail=False):
		vel = self.get_step_velocity(trail)
		return np.sqrt(np.sum(vel**2, axis=1))

	def get_step_data(self,  track_idx=0):
		# return data frame containing information on steps
		_cols = ["track_idx", "sindex", "displacement", "start", "end", "duration"]
		_dtype = [int, int, float, float, float, float]
		n = self.get_nsteps()
		data = {name : np.empty(n, dtype=_dtype[i]) for i, name in enumerate(_cols)}
		data["track_idx"] = np.full(n, track_idx)
		data["sindex"] = np.array(range(n))
		dx = self.get_step_dx()
		data["displacement"] = norm(dx,axis=1)
		_time = self['time'][self.step_idx]
		data["start"] = _time[:-1]
		data["end"] = _time[1:]
		data["duration"] = _time[1:] - _time[:-1]
		data["velocity"]=  data["displacement"] / data["duration"]
		return pd.DataFrame(data)

	def step_angle(self):
		vel = self.get_step_velocity()
		speed = norm(vel,axis=1)
		dv = vel/speed[:,np.newaxis]
		angles = np.arccos(np.sum(dv[1:]*dv[:-1], axis=1))
		return angles

	def decompose_actions_by_velocity(self, tv=0.6, track_idx=0):
		df = self.get_step_data(track_idx=track_idx)
		fast = df["velocity"] > tv
		slow = df["velocity"] < tv
		br = slow[:-1]
		action_idx = np.empty(self.get_nsteps(), dtype=int)
		action_size = np.empty(action_idx.size, dtype=int)
		aidx = 0
		s = 0 
		action_idx[0] = aidx
		last = None
		for i in range(len(br)):
			if br[i]:
				action_size[i-s:i+1] = s+1
				aidx += 1
				s = 0 
				last = i
			else:
				s += 1 
			action_idx[i+1] = aidx
		s = (action_idx.size - last) -1
		action_size[last+1:last+s+1] = s
		df["action_idx"] = action_idx
		df["action_size"] = action_size
		df["large"] = df["action_size"] > 1
		df["fast"] = fast
		df["angle"] = np.insert(self.step_angle(), 0, 0)
		return df

	def decompose_actions(self, arbt=5*np.pi/180, track_idx=0):
		# decompose actions by angle
		df = self.get_step_data(track_idx=track_idx)
		angles = self.step_angle()
		br = angles > arbt
		action_idx = np.empty(len(angles)+1, dtype=int)
		action_size = np.empty(action_idx.size, dtype=int)
		aidx = 0
		s = 0 
		action_idx[0] = aidx
		last = None
		for i in range(len(br)):
			if br[i]:
				action_size[i-s:i+1] = s+1
				aidx += 1
				s = 0 
				last = i
			else:
				s += 1 
			action_idx[i+1] = aidx
		s = (action_idx.size - last) -1
		action_size[last+1:last+s+1] = s
		df["action_idx"] = action_idx
		df["action_size"] = action_size
		df["large"] = df["action_size"] > 1
		return df


	#################################################################################

	# WARN: can't use this procedure by default since it fails on some Fanjin walking tracks
	# ^^ only apply it simulated trajectories then
		
	# IMPORTANT: 2022/04 ~ discovering that the length column and head->trail distance is inconsistent
	# ~ 

	def adjust_projected_axis_by_radius(self, extend):
		## modify x,y,z,trail_x, trail_y,trail_z
		if "width" in self:
			R = np.mean(self["width"])/2
		else:
			R = 0.5
		print("adjusting projected axis using R = {}, extend={}".format(R, extend))
		new = self.copy()
		head = new.get_head()
		trail = new.get_trail()
		center = (head + trail)/2
		segment = head - trail
		xysegment = segment[:,:2]
		xylength = np.linalg.norm(xysegment,axis=1)
		# the factor used to extend the segment
		if extend:
			modfactor = (xylength + 2*R)/xylength
		else: # shrink
			modfactor = (xylength - 2*R)/xylength
		# xylength can be zero, then this division is undefined
		# we set zero here, but it is preferrable to skip over the start of the simulation to avoid vertical initial state
		modfactor = np.nan_to_num(modfactor, 0.0)
		modfactor =  modfactor.clip(0)
		modfactor = modfactor[:,np.newaxis]
		head[:,:2] = center[:,:2] + segment[:,:2]/2 * modfactor
		trail[:,:2] = center[:,:2] - segment[:,:2]/2 * modfactor
		new.set_head(head)
		new.set_trail(trail)
		return new

	def extend_projected_axis_by_radius(self):
		return self.adjust_projected_axis_by_radius(extend=True)

	def shrink_projected_axis_by_radius(self):
		return self.adjust_projected_axis_by_radius(extend=False)

	def copy(self):
		new = TrackLike(self._track.copy(), source=self.source)
		return new


	#################################################################################

	# no slicing
	def get_aspect(self):
		return self['length']/self['width']

	def get_head(self):
		return np.stack([self['x'], self['y'], self['z']], axis=1)

	def get_head2d(self):
		return np.stack([self['x'], self['y']], axis=1)

	def get_trail(self):
		return np.stack([self['trail_x'], self['trail_y'], self['trail_z']],axis=1)

	def get_trail2d(self):
		return np.stack([self['trail_x'], self['trail_y']], axis=1)


	def set_head(self, head):
		self['x'] = head[:,0]
		self['y'] = head[:,1]
		self['z'] = head[:,2]

	def set_trail(self, trail):
		self['trail_x'] = trail[:,0]
		self['trail_y'] = trail[:,1]
		self['trail_z'] = trail[:,2]

	def has_center(self):
		fields = self._track.dtype.fields
		return 'center_x' in fields and 'center_y' in fields

	def get_head_v(self):
		# head velocity
		xy = np.column_stack([self['x'], self['y']])
		dt = self['time'][1:] - self['time'][:-1]
		return (xy[1:] - xy[:-1]) / dt[:,np.newaxis]

	def get_head_d(self):
		xy = np.column_stack([self['x'], self['y']])
		return (xy[1:] - xy[:-1]) 

	def get_center_v(self):
		cxy = self.get_cxy()
		dt = self['time'][1:] - self['time'][:-1]
		return (cxy[1:] - cxy[:-1])/dt[:,np.newaxis]

	def get_cxy(self):
		# centre xy
		if self.has_center():
			return np.stack([self['center_x'], self['center_y']], axis=1)
		else:
			x, y = self['x'], self['y']
			tx, ty = self['trail_x'], self['trail_y']
			cx = (x + tx)/2.
			cy = (y + ty)/2.
			cxy = np.stack([cx, cy], axis=1)
		return cxy

	def get_origin(self):
		if self.has_center():
			return self['center_x'][0], self['center_y'][0], 0.
		else:
			x_ = (self['x'][0]  + self['trail_x'][0])/2.
			y_ = (self['y'][0]  + self['trail_y'][0])/2.
			z_ = (self['z'][0]  + self['trail_z'][0])/2.
			return x_, y_, z_

	def get_frame(self):
		cxy = (self.get_head() + self.get_trail())/2.
		axis = np.stack([self['ax_x'], self['ax_y'], self['ax_z']], axis=1)
		e1 = np.stack([self['rax_x'], self['rax_y'], self['rax_z']], axis=1)
		return cxy, axis, e1
		
	def get_bbox(self):
		x = np.concatenate([self['x'], self['trail_x']])
		y = np.concatenate([self['y'], self['trail_y']])
		return (x.min(), x.max()), (y.min(), y.max())

	def get_bbox_size(self):
		lx, ly = self.get_bbox()
		return np.array([lx[1]-lx[0], ly[1]-ly[0]])


# for experimental data
class TrackLike(Track):

	def __init__(self, track, source='unknown'):
		self._track = track
		self.source = source
		self._init()
		self.step_idx = None


# for convenience in finding track boundaries, etc.

def trxmin(trs, col):
	return min(trmin(trs,col), trmin(trs,'trail_'+col))-0.5
def trxmax(trs, col):
	return max(trmax(trs,col), trmax(trs,'trail_'+col))+0.5

def trmin(trs, col):
	return min([np.min(tr.track[col]) for tr in trs])
def trmax(trs, col):
	return max([np.max(tr.track[col]) for tr in trs])

def maxsize(trs):
	return max((tr.size for tr in trs))
def maxtime(trs):
	return max(tr.track['time'][-1] for tr in trs)
def nptracks(trs):
	return [tr.track for tr in trs]


