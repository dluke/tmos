
import numpy as np
import matplotlib.pyplot as plt

import ctfp3d as tfp
import parameters
import runtfp_global

args = runtfp_global.setup()
#args.Pili.tau = 1.0
cell = tfp.setup_cell(args)
pilus = cell.pili[0]
fstall = args.pget('f_stall')
axvstyle = {'c':'k', 'alpha':0.6, 'linewidth':1., 'linestyle':'--'}

fspace = np.linspace(0, 140, 2000)
def k_dt(f):
    return pilus.calc_k_dt(f)

def tau_eovers(f):
    return pilus.tau_eovers(f)

def plot_k_dt(ax):
    ax.plot(fspace, list(map(k_dt, fspace)))
    ax.plot(0, 0, marker='o')
    ax.plot(args.pget('f_release'), 1./args.pget('tau'), marker='o')
    ax.set_title('rate of detachment vs. force')
    ax.set_xlabel('force (pN)')
    ax.set_ylabel('rate (s^-1)')
    ax.set_yticks([0,10])
    ax.set_yticklabels(['0',r'1/\tau'])
    ax.axvline(fstall,**axvstyle)

def plot_tau_eovers(ax):
    ax.plot(fspace, list(map(tau_eovers, fspace)))
    #ax.set_title(r'\tau_e/\tau_s')
    ax.set_title(r'ratio of extension/retraction for bound pili')
    ax.set_xlabel('force (pN)')
    ax.set_ylabel('rate (s^-1)')
    ax.set_ylim(0,None)
    lfp = args.pget('low_force_proportion')
    hfp = args.pget('high_force_proportion')
    ax.set_yticks([0.,lfp,hfp,1.])
    ax.set_yticklabels(['0.', 'low_force_ratio','high_force_ratio','1.'])
    ax.axvline(fstall,**axvstyle)

fig, axes = plt.subplots(1,2, figsize=(10,6))
ax1, ax2 = axes
plot_k_dt(ax1)
plot_tau_eovers(ax2)
plt.tight_layout()
plt.savefig('force_dependence.png')
plt.show()



