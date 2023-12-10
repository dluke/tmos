import os, sys

single_track_file = 'trackxy.dat'
datdir = 'data'
vtkdir = 'vtk'
plotdir = 'plots'

vtkformat = 'vtk_{:05d}.vtp'
vtkbodyform = 'vtkbody_{:06d}.vtp'
vtkpiliform = 'vtkpili_{:06d}.vtp'



def force_mkdir(dname, verbose=True):
    # avoid triggering an exception 
    if not os.path.exists(dname):
        if verbose:
            print('Directory {} does not exist. Creating it ...'.format(dname))
        os.makedirs(dname)

def setup_vtkdir():
    if not os.path.exists(vtkdir):
        os.makedirs(vtkdir)

# dep eventually
vtkstep = 0
def next_vtk():
    global vtkstep
    vtkf = vtkformat.format(vtkstep)
    vtkstep += 1
    return os.path.join(vtkdir, vtkf)

# vtk output file generator with format as argument
def vtkitf(ddir, form):
    vtki = 0;
    while True:
        yield os.path.join(ddir, form.format(vtki))
        vtki += 1

# special cases
def vtkit(ddir):
    return vtkitf(ddir, vtkformat)

def vtkit_body(): return vtkitf(vtkdir, vtkbodyform)
def vtkit_pili(): return vtkitf(vtkdir, vtkpiliform)

# I need a better way to organise, create and destroy the directories 
# and files we are using to store data
# Can use and object keep track of default locations for storing data
#  and switch between different data organisation methods
#  especially now that we have the capability to produce many .vtp files
class Datapaths(object):
    pass



    #if not os.path.exists(datdir):
        #os.mkdir(datdir)
    ## track output
    #name_out = os.path.join(datdir, 'bacterium_{:05d}.dat')
    ## pili binding/unbinding
    #fpili = os.path.join(datdir, 'pili_{:05d}.dat')


