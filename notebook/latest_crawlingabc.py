# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# copy from crawlingabc.py for latest draft

# %% 
import sys, os
import copy
join = lambda *x: os.path.abspath(os.path.join(*x))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import scipy.stats

from pili import support

import twutils
import pilush
import command
import pili
import parameters
import _fj
import fjanalysis
import twanalyse
import rtw
import abcimplement
from abcimplement import rejection_abc


import pili.publication as pub
import thesis.publication as thesis


# %%
plotting = False

# %%
notedir = os.getcwd()
notename = 'crawlingabc'
root = pili.root
# candidate to compare against
plt.rcParams.update({
	'text.usetex': False,
	'figure.figsize': (20,20),
	'axes.labelsize': 20
	})

# %%
all_idx, ltrs = _fj.slicehelper.load_linearized_trs("all")

# %%
# reference_idx = _fj.load_subset_idx()
cache_dir = "plos_classification/"
ref_idx_list = np.loadtxt(join(cache_dir, "new_crawling.npy")).astype(int)
reference_idx = {'top': ref_idx_list}

objectives = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate']
refdf = fjanalysis.compute_reference_data(ltrs, reference_idx, objectives)
reference = refdf.iloc[0]

# %%
reference

# %%
# updated data is at
# /home/dan/usb_twitching/run/825bd8f/cluster/mc4d
new4d = {}
# new4d["simdir"] = "/home/dan/usb_twitching/run/825bd8f/cluster/mc4d"
new4d["simdir"] = "/home/dan/usb_twitching/run/latest/mc4d"
new4d["objectives"] = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate', 'fanjin.top.ks_statistic', 'nbound.mean', 'ntaut.mean']
new4d = abcimplement.load_problem_simulation(new4d)

# %%
bound_pili_participation = new4d["data"].get("bound_pili_participation")
# bound_pili_participation[]
bound_pili_participation[_accepted.index].mean(), bound_pili_participation.mean()

# %%
# print problem
nsamples = int(1e4)
N = 50
print("accept {}/{}".format(N,nsamples))

# %%
# one statistic at a time
new4d["params"] = new4d["data"].paramsdf(new4d["objectives"])
# abcimplement.transform_anchor_parameter_data(new4d)
abcimplement.transform_vret_parameter_data(new4d)
_objectives = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate']
statdf, statref = abcimplement.regularise_stats(new4d["params"], reference, _objectives)
print('marginal standard deviations')
statdf.attrs['std']

# %%
if plotting:
	for objective in new4d["objectives"]:
		_regdf = statdf[new4d["problem"]["names"] + [objective]]
		_accepted = rejection_abc(_regdf, [objective], statref, N)
		# rename = {k:k for k in _accepted.keys()}
		# rename["anchor_angle_smoothing_fraction"] = "anchor
		# _accepted.rename(columns=rename, inplace=True)
		abcimplement.problemplot4d(new4d["problem"], _accepted, objective)
		plt.tight_layout()

# %%
# combine lvel, deviation.var, activity

# all three remaining simple metrics
# _objectives = ["lvel.mean", "deviation.var", "fanjin.top.ks_statistic"]

# _objectives = ["lvel.mean", "deviation.var", "ahat.estimate"]

_objectives = ["lvel.mean", "deviation.var", "qhat.estimate", "ahat.estimate"]

_regdf = statdf[["uid"] + new4d["problem"]["names"] + _objectives]
_accepted = rejection_abc(_regdf, _objectives, statref, N)
m_par, v_par = abcimplement.mean_accepted(new4d["problem"], _accepted)
abcimplement.describe_abc(new4d, _accepted)
prime_accepted = _accepted

df = abcimplement.tabulate_inference(new4d["problem"], _accepted, _objectives)
m_par = {k : v for k, v  in zip(df["parameter"], df["MLE"])}

df


# %%
_accepted.corr()

# %%
params = new4d["params"]
params

# %%

accparams = new4d["params"].iloc[_accepted.index]
accparams["score"] = _accepted["score"]


# %%
# compute the confidence interval for nbound, ntaut
accparams
accparams["nbound.mean"].min(), accparams["nbound.mean"].max()

epsilon = _accepted['score'].iloc[-1]
scores = _accepted['score'].to_numpy()
kernel = abcimplement.new_epanechnikov(epsilon)
weights = kernel(scores)
nbins = 10

fig, ax = plt.subplots(figsize=(5,5))
sns.histplot(x=accparams['ntaut.mean'], bins=nbins, ax=ax, color='#88CB89', weights=weights, 
		alpha=0.6, kde=False,  line_kws={'linewidth':4})

fig, ax = plt.subplots(figsize=(5,5))
sns.histplot(x=accparams['nbound.mean'], bins=nbins, ax=ax, color='#88CB89', weights=weights, 
		alpha=0.6, kde=False,  line_kws={'linewidth':4})

# %%
lpar = [r'$\tau_{\mathrm{dwell}}$', r'$\kappa$', r'$v_{\mathrm{ret}}$', r'$k_{\mathrm{spawn}}$']
fig, axes = abcimplement.perfectplot4d(new4d["problem"], _accepted, mpar=m_par, lpar=lpar)
pub.save_figure('top_rejectionabc', notename)

# %%

df = abcimplement.tabulate_inference(new4d["problem"], _accepted, _objectives)
pub.save_dfhtml(df, "fj_crawlingabc", notename)
df
# list(zip(df['parameter'], df['MLE']))

# %%
# pull the best simulation
import readtrack

# udir = join(new4d["simdir"], new4d["lookup"][0][_accepted.index[0]])
udir = join(new4d["simdir"], new4d["lookup"][0][_accepted.index[1]])
udir

# %%
_accepted

# %%
# twutils.sync_directory(udir)
args = parameters.thisread(directory=udir)
_trs = readtrack.trackset(ddir=join(udir, 'data/'))
bestsim = [_fj.linearize(tr) for tr in _trs]

# %%

# best parameters
ref_ltrs = [ltrs[idx] for idx in ref_idx_list]
ref_vel = np.concatenate([ltr.get_step_speed() for ltr in ref_ltrs])
_accepted.iloc[:10]
accparams.iloc[:10]

# %%

with plt.rc_context({'text.usetex': True, 'font.size':28}):
	fig, ax = plt.subplots(figsize=(5.5,5))
	sim_vel = np.concatenate([ltr.get_step_speed() for ltr in bestsim])

	# xlim = (0, np.quantile(sim_vel, 0.98))
	xlim = (0, 1.8)
	xlim = (0, 1.0)

	ks, p = scipy.stats.ks_2samp(ref_vel, sim_vel)
	print("ks statistic = {:0.3f}".format(ks))
	hstyle = {'stat':'density', 'alpha':0.3, 'fill':True, 'element':'step'}
	common = {"bins": 50, 'binrange':xlim}
	common.update(hstyle)
	sns.histplot(ref_vel, ax=ax, label="reference", **common)
	sns.histplot(sim_vel, ax=ax, color="orangered", label="simulation", **common)
	ax.set_xlim(xlim)
	ax.legend(fontsize=28)
	ax.set_xlabel(r"step velocity (\textmu m/s)", fontsize=30)
	ax.set_ylabel("density", fontsize=30)
	ax.grid(False)   

plt.tight_layout()
pub.save_figure("best_sim_crawling_lvel_distrib", notename)


# %%
udir
# _u_DbaHFe2d

# %%
# stacked bars
with command.chdir(join(udir, 'localrun')):
	ptdata = readtrack.piliset()[0]


# %%
visible_ranges = [(0,99.0), (0.3,99.0),(1.0,99.0),(2.0,99.0)]

timelime, npy, npy_in_range, bound_in_range, taut_in_range = pilush.count_visible(ptdata, visible_ranges=visible_ranges)
clrange = [npy_in_range, bound_in_range, taut_in_range]


mean_npili = np.mean(npy)
form = 'Average number of pili in range {} is {:3.1f}/{:3.1f}'
for vr, nlist in list(npy_in_range.items()):
    print((form.format(vr, np.mean(nlist), mean_npili)))
print()

mean_npili = np.mean(bound_in_range[visible_ranges[0]])
for vr, nlist in list(bound_in_range.items()):
    print((form.format(vr, np.mean(nlist), mean_npili)))
print()


mean_npili = np.mean(taut_in_range[visible_ranges[0]])
for vr, nlist in list(taut_in_range.items()):
    print((form.format(vr, np.mean(nlist), mean_npili)))




# %%

with mpl.rc_context({'font.size': 20, 'text.usetex': True}):
	fig, ax = plt.subplots(figsize=(3,4))
	# support.plot_stacked_bars(ax, clrange, visible_ranges)
	
	support.plot_bars(ax, clrange, visible_ranges)

pub.save_figure("crawling_best_accepted_tfp_count_bars.png", notename)
thesis.save_figure("crawling_best_accepted_tfp_count_bars.png", notename)


# %%
df

# %%
row = df.loc[df["parameter"] == "pilivar"]
triplet = [row["confidence(0.05)"].iloc[0], row["MLE"].iloc[0], row["confidence(0.95)"].iloc[0]]
fig, ax = support.new_plot_tfp_distrib(triplet)
pub.save_figure("crawling_predicted_kappa", notename)

# %%
# load the xydistrib from the best simulation
xydisp = np.load(join(udir, "xydisp.npy"))
xydisp = xydisp[xydisp > 0]

hstyle = {'stat':'density', 'alpha':0.3, 'fill':True, 'element':'step'}
with plt.rc_context({'text.usetex': True, 'font.size':20}):
	fig, ax = plt.subplots(figsize=(5.5, 5))
	def plot_pdisp_distrib(ax, xydisp, hstyle):
		sns.histplot(xydisp, ax=ax, **hstyle)
		vstyle = dict(c='k', alpha=0.6,  linestyle='--')
		ax.axvline(0.12, **vstyle)
		ax.set_xlabel(r"per TFP displacement $(\mu m)$")
		ax.set_xlim((0.,0.3))
	plot_pdisp_distrib(ax, xydisp, hstyle)
best_xyd = np.mean(xydisp)

# %%

simdir = new4d["simdir"]
def save_accepted_samples(simdir):
	accept_uid = _accepted["uid"].to_numpy()
	cachefile = join(simdir, "standard_accepted_uid.npy")
	print("saving accepted uid list to ", cachefile)
	np.savetxt(cachefile, accept_uid, fmt="%s")
save_accepted_samples(simdir)


# %%

# pull the xydisp.npy data for the accepted samples
def load_accepted_uid(simdir):
	return np.loadtxt(join(simdir, "standard_accepted_uid.npy"), dtype=object)

def load_accepted_xydisp(simdir):
	uidlist = load_accepted_uid(simdir)
	data = [np.load(join(simdir, uid, 'xydisp.npy')) for uid in uidlist]
	data = list(map(lambda arr: arr[arr>0], data))
	return data

# load_accepted_uid(simdir)
uidlist = load_accepted_uid(simdir)
xydispdata = load_accepted_xydisp(simdir)

# %%
# target = join(pili.root, '../run/825bd8f/target/t0_anchor')
# TODO
target = join(pili.root, "../run/825bd8f/line_search/alpha_search")
dc = rtw.DataCube(target=target)
xydata = [np.load(join(directory, "xydisp.npy")) for directory in dc.dircube]
pdisp = [np.mean(xyd[xyd> 0]) for xyd in xydata]


# %%
means = [np.mean(xyd) for xyd in xydispdata]
np.mean(means)

# %%
x = _accepted["kb_sh"].to_numpy()
with mpl.rc_context({'text.usetex':True, 'font.size':24}):
	fig, ax = plt.subplots(figsize=(4,4))
	mscat = {'edgecolor':'white', 'linewidth' : 1.0, 'alpha' : 0.8, 's': 60}
	ax.scatter(x, means, **mscat)
	ax.set_ylim((0,None))
	ax.set_xlim((0,1))
	# ax.set_xticks([0,np.pi/4,np.pi/2])
	# ax.set_xticklabels(['0',r'$\pi/4$',r'$\pi/2$'])

	# ax.scatter(m_par["anchor_angle_smoothing_fraction"], best_xyd)

	m0 = {'c':'indianred', 'marker':'D', 'alpha':1.0, 'markersize':4}
	base = np.array(dc.basis[0]) * np.pi/2
	ax.yaxis.set_major_locator(plt.MaxNLocator(5))

	# TODO
	# ax.plot(base, pdisp, **m0)

	# ax.scatter(_accepted.iloc[0]["anchor_angle_smoothing_fraction"], best_xyd, marker='X', edgecolor='white', s=140)

	ax.legend(["accepted", "line search"], loc='upper right', bbox_to_anchor=(1.18, 0.99), fontsize=18, framealpha=1.0)

	ax.set_ylabel(r'$\langle \Delta x \rangle$ (\textmu m)')
	vret = r'$v_{\mathrm{ret}}$'
	ax.set_xlabel(vret)

	inset = ax.inset_axes([0.52, 0.24, 0.6, 0.4])
	def _plot(ax, xydisp, hstyle):
		vstyle = dict(c='k', alpha=0.6,  linestyle='--')
		orange = "#FF7F0E"
		sns.histplot(xydisp, ax=ax, color=orange, **hstyle)
		ax.axvline(0.12, **vstyle)
		ax.set_xlabel(r"$ \Delta x$ (\textmu m)", fontsize=16)
		ax.set_ylabel(r"density", fontsize=20)
		ax.set_xlim((0.,0.3))
		# ax.annotate(r'$\delta_{\text{step}}', (0.4,0.5), xytext=(0.6, 0.5), textcoords='axes fraction')
		ax.annotate(r'$\delta_{\mathrm{step}}$', (0.12, 10), xytext=(0.6, 0.5), textcoords='axes fraction', fontsize=16, arrowprops=dict(facecolor='black', width=1, shrink=0.05, headwidth=4))

		# ax.annotate('{:.3f}'.format(ratio), (ratio,0), xytext=(0.7,0.5), textcoords='axes fraction', arrowprops=dict(facecolor='black', width=2, shrink=0.05)) 

	_plot(inset, xydisp, hstyle)

pub.save_figure('crawling_accepted_mean_pdisp', notename)
thesis.save_figure('crawling_accepted_mean_pdisp', notename)

# %%
# we are going to create a simulation at this value and then vary \alpha
print("alpha", m_par["anchor_angle_smoothing_fraction"] / (np.pi/2) )
m_par

# %%

anchor, anchor_contraction = _get_par(dc_anchor, dlstat)
anchor, a_disp = _get_par(dc_anchor, dispstat)
anchor *= np.pi/2
print('anchor', anchor)
#
