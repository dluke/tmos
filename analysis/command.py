import os
import sys
import types
import collections
import time

import datapaths

def options(local):
    opt = {}
    for fu, loc in list(local.items()):
        if isinstance(loc, types.FunctionType) and not hasattr(loc, 'exempt'):
            opt[fu] = loc
    return opt


def process(local):
    args = sys.argv[1:]
    nargs = len(args)

    local = options(local)
    # No function given so print options and exit
    if len(args) == 0:
        for fu in sorted(local.keys()):
            print(fu)
        raise RuntimeError('No function specified. Exiting.')

    fname = args[0]
    fargs = args[1:]
    try:
        f_using = local[fname]
    except KeyError:
        print('I don\'t have a function \'%s\'', fname)
        raise

    if len(fargs) == 0:
        #print 'No arguments given to ' , fname

        if hasattr(f_using, 'defaults'):
            fargs = f_using.defaults
            print('Using defaults ', f_using.defaults)
        else:
            fargs = []

    ff = fargs
    for i, f in enumerate(ff):
        try:
            ff[i] = eval(f)
        except:
            # This string has a '.' and look like an object
            #ff[i] = f
            pass

    return f_using, ff

# IO utility

#print square data to file, first column is int and rest are floats.
def dump(dd, fo, htag='#', outstr=None):
    nc = len(list(dd.keys()))
    #fo.write(''.join([htag+' ', '%s\t '*nc, '\n']) % tuple(dd.keys()))
    fo.write((htag+' '+'{:<14} '*nc + '\n').format(*list(dd.keys())))
    ddv = list(dd.values())
    nr = len(ddv[0]) # assumption
    if not outstr:
        outstr = '      '+'{:<14} '*nc + '\n' 
    for i in range(nr):
        dv = [d[i] for d in ddv]
        fo.write(outstr.format(*dv))

# no 'read_dump' method?

def readdump(f):
    table = collections.OrderedDict()
    pass 

# convenience methods for writing meta data int the same directory as the plots

def metafile(form):
    return os.path.splitext(form)[0]+'.meta'

lineform = '{}\t\t{}\n'
def writemeta(meta, forfile):
    # meta is a dictionary of metadata
    mout = metafile(forfile)
    print('writing metadata to ', mout)
    with open(mout, 'w') as fm:
        for key, v in list(meta.items()):
            fm.write(lineform.format(key, v))


## From SAMoS python. descriptor for organising analysis plots and data.

import os
from matplotlib import pyplot as plt

def mkdir(dirpath):
    os.system('mkdir -p {}'.format(dirpath))

# helper for defaultsave
pdir = 'plots/'
ddir = os.path.join(pdir, 'data/')
def result_dirs(pdir, ddir):
    if not os.path.exists(pdir):
        print('make directory {}'.format(pdir))
        mkdir(pdir)

    if not os.path.exists(ddir):
        mkdir(ddir)

def default_result_dirs():
    result_dirs(pdir, ddir)

def save_all_figures(f, pdir=pdir, extra='', svg=False, noback=False):
    print('save all figures...')
    for f_id in plt.get_fignums():
        figure = plt.figure(f_id)
        base = f.__name__
        if len(extra) != 0:
            base += '_' + extra
        if hasattr(figure, '_plot_name'):
            base+= '_'+figure._plot_name
        pbase = os.path.join(pdir, base)
        out = pbase+'.svg'
        outpng = pbase+'.png'
        print('writing plot to ', outpng)
        plt.savefig(outpng, format='png', transparent=noback)
        if svg:
            plt.savefig(out, format='svg', transparent=noback)


from functools import wraps
# descriptor for saving the plots to a default name under plots/ and the data under plots/data
def defaultsave(autoshow=False, pdir=pdir, ddir=None, extra='', svg=False, noback=False):
    if ddir == None:
        ddir = os.path.join(pdir, 'data/')
    def rdec(f):
        @wraps(f)
        def saved(*args, **kw):
            ret = f(*args, **kw)
            outd = ret

            # make the target directories if they don't exist
            result_dirs(pdir, ddir)
            save_all_figures(f, pdir, extra, svg, noback)

            # if outd:
                # datout = os.path.join(ddir, base+'.dat')
                # print('saving data to ', datout)
                # with open(datout, 'w') as fo:
                #     dump(outd, fo)

            if autoshow:
                plt.show()
            return ret
            #plt.show()
        return saved
    return rdec

# little helper for saving plot while sending notification to stdout.  
def saveplt(base, pdir=pdir, svg=False, transparent=False):
    datapaths.force_mkdir(pdir)
    pbase = os.path.join(pdir, base)
    print(('writing plot to ', pbase+'.png'))
    plt.savefig(pbase+'.png', format='png', transparent=transparent)
    if svg:
        plt.savefig(pbase+'.svg', format='svg', transparent=transparent)
    
def savefig(fig, base, pdir=pdir, svg=False, transparent=False):
    datapaths.force_mkdir(pdir)
    if base.endswith('.png'):
        base = base[:-4]
    pbase = os.path.join(pdir, base)
    print(('writing plot to ', pbase+'.png'))
    fig.savefig(pbase+'.png', format='png', transparent=transparent)
    if svg:
        fig.savefig(pbase+'.svg', format='svg', transparent=transparent)

"""
take an iterator that returns arguments for the wrapped function
along with an index. 
Call the function with the argument and then plt.savefig to
a unique file.
imap is a function that maps the plot index i to the track index

TODO replate f(ft) with iterator and pass a plotting function in advance
"""
def defaultsaveall(npyiter, imap=lambda x: x, pdir=pdir, noback=False, nform='{:05d}'):
    def rdec(f):
        formdir = os.path.join(pdir, f.__name__)
        datapaths.force_mkdir(formdir)
        form = os.path.join(formdir, '_'.join(['line', nform+'.png']))
        metaf = os.path.join(formdir, 'meta_map.txt')
        if os.path.exists(metaf):
            os.remove(metaf)
        @wraps(f)
        def saveall():
            for i,ft in npyiter:
                plt.clf()
                f(ft)
                out = form.format(i)
                print('writing to', out)
                plt.savefig(out)
                with open(metaf, 'a') as mf:
                    mf.write(nform+' '+nform+'\n'.format(i, map(i)))
            return 
        return saveall
    return rdec

# for running code in a particular directory
class chdir(object):
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self._curr_path = os.getcwd()
        os.chdir(self.path)

    def __exit__(self, type, value, tb):
        os.chdir(self._curr_path)

# for more descriptive jupyter notebooks
class loadexplain(object):
    def __init__(self, message):
        self.message = message

    def __enter__(self):
        print('loading {}...'.format(self.message))

    def __exit__(self, type, value, tb):
        print('finished loading {}.'.format(self.message))


class timer(object):

    def __init__(self, message, fmt="{:.3f}"):
        self.message = message
        self.fmt = fmt

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, tb):
        end = time.time()
        diff = end - self.start
        form = "{} (%ss)" % self.fmt
        print(form.format(self.message, diff))
