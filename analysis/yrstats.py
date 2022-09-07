

import os, itertools
import numpy as np
import matplotlib.pyplot as plt

import command, filesystem

import readmat
import matdef
import fjstats
import readyr
import _analysefj

pdir = readyr.plots
partial = readyr.partial

yrdenoiseddir = '/home/dan/usb_twitching/yr/denoised_1'
yrdisp = os.path.join(yrdenoiseddir, 'dr_0.1s_denoised_all')

# 
def load_vel():
    print('loading displacements from ', yrdisp)
    disp = np.loadtxt(yrdisp)
    # cut to 0.5 \mu m maximum
    vel = disp*10
    vel = fjstats.cut_vel(vel)
    return vel


@command.defaultsave(pdir=pdir)
def view_vel():
    vel = load_vel()
    outd = plotutils.kdeplot(vel)
    print('average velocity', np.mean(vel))
    #plt.hist(disp, bins=100) 
    plt.ylabel('probability distribution')
    plt.xlabel('0.1 second velocity')
    plt.tight_layout()
    return outd


# READING
def _read_yt():
    # read high time res data for analysis with this module
    yt = readyr.inithres()
    # seriously, what is up with YR angle?
    yt.override_orientation_using_lead_trail()
    return yt

   
# debugging test columns are read correctly
def print_first_line():
    #yt = readyr.inithres()
    # inithres
    highresfile = '/home/dan/twitching/yr/190314/tr.txt'
    pix_to_mu = 0.042
    # no conversions
    yt = readmat.ExpTracks.useyr(highresfile, timestep=0.1,
            colmap=matdef.yrcolmap2, yrpix_to_mu=pix_to_mu, convert=False)
    track = yt.tracks[0]
    # print a bunch of stuff, mainly the first 2 values of each column, columns ordered
    print(yt.colmap)
    rcolmap = {v:k for k,v in list(yt.colmap.items())}
    print(track.shape)
    for value  in range(len(yt.colmap)):
        col = rcolmap[value]
        trcol = yt.get_col(track, col)
        print(value, col, trcol[0], trcol[1])
    # more debugging
    orient = yt.get_whole_col('orientation')
    #print 'orientation', (180./np.pi) * np.array([orient.min(), orient.max()])
    print('orientation',  np.array([orient.min(), orient.max()]))


def check_against_yr():
    # can also check by eye now that I corrected my mistake with reading columns
    yt = readyr.inithres()

    # seriously, what is up with YR angle?
    yt.override_orientation_using_lead_trail()
    # get the fast percent for speed > 0.01 \mu m

    col = matdef.LEAD
    pu = _analysefj.disps(yt, timestep=yt.timestep, cols=col)/yt.timestep
    fastpu= pu[pu>0.2]
    # check absolute value of velocity against yr
    print(fastpu.shape)
    print(np.min(pu), np.max(pu))
    #fast = np.percentile(pu, fastpercent)

###################################################################################
# axis/head velocity deviation angle in YR style
# YR style is 
def frtheta(yt, fastpercent):
    # 
    col = matdef.LEAD
    pu = _analysefj.disps(yt, timestep=yt.timestep, cols=col)/yt.timestep
    fast = np.percentile(pu, fastpercent)

    # compute FanJin reorientation angles
    thetas = np.abs(_analysefj.vel_dot_ax(yt, timestep=yt.timestep, point=col))
    print('average R angle', to_degrees(np.mean(thetas))) # should be zero
    #print thetas.min(), thetas.max()

    fast_th= thetas[pu>=fast]
    print('Plotting angles using {}% of data'.format(100.*fast_th.size/thetas.size))
    print('no. values {}'.format(fast_th.size))
    print('Using actions with velocity > {}'.format(fast))
    deg_thetas = 180./np.pi * fast_th
    return deg_thetas

def fr(yt, fastpercent):
    deg_thetas = frtheta(yt, fastpercent)
    plt.title("Reorientations Fast {}%".format(100-fastpercent))
    histstyle = {'density':True, 'bins':10, 'alpha':0.4}
    plt.hist(deg_thetas, **histstyle)
    yr_angle(deg_thetas)

@command.defaultsave(pdir=pdir)
def yr_angle(thetas):
    plotutils.kdeplot(thetas)
    plt.xlabel("\"reorientation angle\"")
    plt.ylabel("probability density")
    plt.tight_layout()
 
# main method deviation/reorientation angle
def yrdeviation():
    yt = readyr.inithres()
    fastpercent = 99
    # seriously, what is up with YR angle?
    #yt.override_orientation_using_lead_trail()  # now done automatically 
    # get the fast percent for speed > 0.01 \mu m
    fr(yt, fastpercent=fastpercent)


# helpers
def saveas(form, i, sdir='./', ext=''):
    out= os.path.join(sdir, form.format(i)) + ext
    filesystem.safemkdir(out, verbose=True)
    plt.savefig(out)


def to_degrees(a):
    return 180./np.pi * a

def get_delta_bangle(yt):
    orientation = yt.get_cols('orientation') # returns list
    def get_delta(trackor):
        diff = trackor[1:] - trackor[:-1]
        for i in range(diff.size):
            if np.abs(diff[i]) > np.pi: # assume this doesn't happen naturally
                diff[i] = 2*np.pi - np.abs(diff[i])
        return diff
    return [to_degrees(get_delta(x)) for x in orientation]

#################################################################################
# experimental

#axislims = [(-200,200), (0,0.1), (-0.5,0.5)]
axislims = [(-200,200), (0,0.1), (None,None)]
def fjr_corr(yt, form, axislims=axislims):
    """
    yt : Exptracks object
    form : is output form, i.e. yrplots/singletrack/r_angle_{:04d}.png'
    """
    # sequence of per track, R angle, displacement and \delta body angle 
    # 

    figure, axes = plt.subplots(3,1, sharex=True, figsize=(15,10))
    ax, axv, axb = axes
    ylabels = ['R angle\n (degrees)', 'displacement\n (\mu m)', 
            '\delta body angle\n (degrees)']
    colors = itertools.cycle(['b', 'g', 'r'])

    # compute displacements and angles
    col = matdef.LEAD
    disps = _analysefj.track_disps(yt, yt.timestep, cols=col)
    pertrack = _analysefj.track_vel_dot_ax(yt, yt.timestep, point=col)
    # change in body axis angle
    delta_bangle = get_delta_bangle(yt)

    for axis, ylims in zip(axes, axislims):
        axis.set_ylim(*ylims)
    for i, track in list(enumerate(yt.tracks)):
        for axis in axes:
            axis.clear()
        ddotax = pertrack[i]
        ddotax = to_degrees(ddotax)
        # convert to degrees
        #style = {'marker':'o', 'linestyle':'--', 'linewidth':1, 'markersize':4}
        style = {'linestyle':'--', 'linewidth':1, 'markersize':4}
        #
        ax.plot(ddotax, color=next(colors), **style)
        #
        axv.plot(disps[i], color=next(colors), **style)
        #
        axb.plot(delta_bangle[i], color=next(colors), **style)
        #axb.set_ylim(-3.0, 3.0)

        #
        for j, axis in enumerate(axes):
            axis.set_ylabel(ylabels[j])
        axb.set_xlabel('frame (interval 0.1 s)')

        plt.tight_layout()
        saveas(form, i, sdir='.')

def ordered_corr(yt):
    # plotting to visualise the correlation between R angle and displacement
    # 
    col = matdef.LEAD
    pu = _analysefj.disps(yt, timestep=yt.timestep, cols=col)/yt.timestep

    # compute FanJin reorientation angles
    thetas = to_degrees(np.abs(_analysefj.vel_dot_ax(yt, timestep=yt.timestep, point=col)))
    sortidx = np.argsort(thetas)
    spu = np.clip(pu[sortidx], 0., 1.)
    sthetas = thetas[sortidx]

    style = {'marker':'o', 'linestyle':'None', 'markersize':2}
    plt.plot(spu, sthetas, **style)
    plt.show()

# naive plot the distribution of angle * vel
@command.defaultsave(pdir=pdir)
def plot_fjr_corr(distrib):
    #distrib = distrib[distrib<0.1] # failed to identify high velocity 
    plotutils.kdeplot(distrib)
    plt.xlabel("vel * angle")
    plt.ylabel("probability distribution")
    plt.tight_layout()

def yr_steptime(yt, step_d=0.24, scaling=1., fit=False, linekw={}):
    #step_d = 0.042 # one pixel ...
    step_d = step_d
    fjstats.fj_steptime(yt, step_d, scaling, fit, linekw)
    plt.savefig(os.path.join(pdir, '_stime_distribution.png'))

@command.defaultsave(pdir=pdir)
def divider_dimension(ft):
    #dbasis = np.linspace(0.03, 0.60, 21, True)
    #dbasis = np.linspace(0.03, 1.20, 41, True)
    dbasis = np.logspace(-2, 2, 41, True)
    fjstats._divider_dimension(ft, dbasis)

### apply yrstats functions to fj

import _fj  # fj reader
def apply_fjr_corr():
    # would like to see the same timescale
    ft = _fj.npyloadall(list(range(100)))
    form = 'fjplots/singletrack/r_angle_{:04d}.png'
    fjr_corr(ft, form)

if __name__=='__main__':
    # debugging
    #view_vel()
    #check_against_yr()
    
    yt = _read_yt()

    # METRICS
    # R angle metric
    yrdeviation()
    #yr_steptime(yt)

    # Visualisation
    #  
    #form = 'yrplots/singletrack/r_angle_{:04d}.png'
    #fjr_corr(yt,form)
    #
    #ordered_corr(yt)


    #divider_dimension(yt)


