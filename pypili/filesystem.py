
import os, sys
from glob import glob

import collections

# process directories
exclude_dirs = ['plots', 'data', 'save', 'vtk', 'notes', 'tmp', 'extra', 'save_plots']

def thisdir():
    return os.path.basename(os.getcwd())

# make the directory of path if it doesn't exist
def safemkdir(path, verbose=False):
    dirpath = os.path.dirname(path)
    try:
        if verbose:
            print("making directory for file {}".format(path))
        os.makedirs(dirpath)
    except OSError:
        if not os.path.isdir(dirpath):
            raise


"""
make this inclusive instead of exclusive
"""
def getdirs(target='.'):
    allparams = read_param_meta(target)
    alldirs = sorted( next(os.walk(target))[1] )
    def is_exclude(d):
        # return any([ex in d.split('/') for ex in exclude_dirs])
        return not any([d.startswith(p_name) for p_name in allparams])
    ddirs = [dr + os.sep for dr in alldirs if not is_exclude(dr)]
    return ddirs

def read_param_meta(target='.'):
    """
    return parameters ordered dict
    """
    falp = os.path.join(target, 'parameters.meta')
    if not os.path.exists(falp):
        raise IOError("Can't find {}".format(falp))
    allparams = collections.OrderedDict()
    with open(falp, 'r') as fp:
        fp.readline()
        fp.readline()
        for line in fp:
            sl = line.split(' ')
            name, s_value = sl[0], ''.join(sl[1:]).strip()
            try:
                parls = eval(s_value)
            except:
                raise ValueError("failed to convert parameter list string to python list\n{}"
                        .format(s_value))
            assert(type(parls) == type([]))
            allparams[name] = parls
    return allparams

def par_split(st, parnames):
    """
    split a simulation folder into parameter names and values
    Using meta data from parameters.meta
    """
    if st.endswith('/'):
        st = st[:-1]
    split = []
    for name in parnames:
        if st[0] == '_':
            st = st[1:]
        l1 = st.find(name)
        if (l1 != 0):
            raise RuntimeError("Start of dirname should be param name.\n{}\n{}"
                    .format(st, name))
        split.append(name)
        st = st[len(name+'_'):]
        v1 = st.find('_')
        if v1 == -1: v1 = None
        split.append(st[:v1])
        st = st[v1:]
    return split

# does this need to exist? I use getdirs everywhere
def procdirs(target='.'):
    ddirs = getdirs(target)
    print('found directories', ddirs)

    # parameter names without '_'
    #pars, val = zip(*[dd.split('_') for dd in ddirs])
    # parameters names with '_' using parameters.meta metadata
    print(list(zip(*[par_split(dd) for dd in ddirs])))
    pars, val = list(zip(*[par_split(dd) for dd in ddirs]))
    assert pars[:-1] == pars[1:]
    return pars[0], ddirs, val

# replaced by scripts/datawalk scipts/seqwalk
def operdirs(anything, *args):
    _, ddirs, _ = fs.procdirs()
    cwd = os.getcwd()
    for i, dd in enumerate(ddirs):
        os.chdir(dd)
        anything(*args)
        os.chdir(cwd)



