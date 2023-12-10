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
# Rejection ABC against a simulated target

# %%
import pili.publication as pub
import twutils
import abcimplement
import sobol
import rtw
import twanalyse
import fjanalysis
import _fj
import parameters
import pili
import stats
import command
import readtrack
import scipy.stats
import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys
import os
import copy
join = lambda *x: os.path.abspath(os.path.join(*x))


# %%
plotting = False

# %%
# config
plt.rcParams.update({
    'text.usetex': False,
    'figure.figsize': (20, 20),
    'axes.labelsize': 16
})
notename = 'simabc'
verbose = False


# %%
# simtarget = "/home/dan/usb_twitching/run/825bd8f/target/t0"
simtarget = "/home/dan/usb_twitching/run/latest/target/t0"
with command.chdir(simtarget):
    ltarget = stats.load()
    args = parameters.thisread()
_simobjective = ['lvel.mean', 'deviation.var', 'qhat.estimate',
                 'ahat.estimate', 'kmsd.mean', 'nbound.mean', 'ntaut.mean']
simref = {name: twutils.make_get(name)(ltarget) for name in _simobjective}
_interest = ['dwell_time', 'k_spawn', 'pilivar', 'kb_sh']

simpar = {par: args.pget(par) for par in _interest}
simpar['kb_sh'] *= 0.004
simpar, simref

# %%

_pdisp = np.load(join(simtarget, "xydisp.npy"))
_pdisp = _pdisp[_pdisp > 0]
simref["pdisp"] = np.mean(_pdisp)


# %%
sim4d = {}
sim4d["simdir"] = "/home/dan/usb_twitching/run/latest/mc4d"
sim4d["objectives"] = ['lvel.mean', 'deviation.var', 'qhat.estimate',
                       'ahat.estimate', 'sim.t0.ks_statistic', 'kmsd.mean', 'nbound.mean', 'ntaut.mean']
sim4d = abcimplement.load_problem_simulation(sim4d)
sim4d["problem"]

# %%
def load_pdisp(simdir):
    data = [np.load(join(simdir, uid, 'xydisp.npy')) for uid in sim4d["lookup"][0]]
    data = list(map(lambda arr: arr[arr>0], data))
    pdisp = list(map(np.mean, data))
    return pdisp
pdisp = load_pdisp(sim4d["simdir"])
assert(len(sim4d["lookup"][0]) == len(pdisp)) #! careful

# %%

# %%
# ABC config
N = 50
print('{}/{}'.format(N, sim4d["M"]))

# %%
# one statistic at a time
_objectives = sim4d["objectives"]
sim4d["params"] = sim4d["data"].paramsdf(_objectives)
sim4d["params"]["pdisp"] = pdisp
abcimplement.transform_vret_parameter_data(sim4d)
statdf, statref = abcimplement.regularise_stats(
    sim4d["params"], simref, _objectives)

# %%
N = 50
# all four simple metrics
warnings.filterwarnings("ignore")

# _objectives = ["sim.t0.ks_statistic", "deviation.var", "qhat.estimate", "ahat.estimate"]
# _objectives = ["lvel.mean", "deviation.var", "nbound.mean", "ahat.estimate"]
# _objectives = ["lvel.mean", "deviation.var", "nbound.mean"]
# _objectives = ["lvel.mean", "deviation.var", "ahat.estimate", "qhat.estimate", "pdisp"]
# _objectives = ["lvel.mean", "pdisp"]
# _objectives = ["lvel.mean", "deviation.var", "qhat.estimate", "ahat.estimate", "ntaut.mean"]
# _objectives = ["lvel.mean", "deviation.var", "ahat.estimate", "qhat.estimate", "ntaut.mean"]
_objectives = ["lvel.mean", "deviation.var", "qhat.estimate", "ahat.estimate"]

_regdf = statdf[['uid'] + sim4d["problem"]["names"] + _objectives]
_accepted = abcimplement.rejection_abc(_regdf, _objectives, statref, N)

df = abcimplement.tabulate_inference(
    sim4d["problem"], _accepted, _objectives, simpar=simpar)
m_par = {k: v for k, v in zip(df["parameter"], df["MLE"])}

prime_accepted = _accepted
pub.save_dfhtml(df, "sim_crawlingabc", notename)
df

# %%
vret = r'$v_{\mathrm{ret}}$'
lpar = [r'$\tau_{\mathrm{dwell}}$', r'$\kappa$', r'$v_{\mathrm{ret}}$', r'$k_{\mathrm{spawn}}$']
if True:
    fig, axes = abcimplement.perfectplot4d(
        sim4d["problem"], _accepted, simpar=simpar, mpar=m_par, lpar=lpar)
    pub.save_figure('new_sim_crawling_abc', notename, fig)


# %%
_accepted = abcimplement.rejection_abc(_regdf, _objectives, statref, N)

names = sim4d["problem"]["names"]
_accepted[names].corr()

# %%
# TODO if pdisp is relatively constant while ntaut/kspawn/etc vary, 
# it means that per pilus displacement is driving the observable statistics to be similar.
params = sim4d["params"]
accparams = params.iloc[_accepted.index]
accparams["score"] = _accepted["score"]
# accparams

# %%
params[names + ['pdisp']].corr()

# %%
accparams["pdisp"].mean()

# %%

fig, ax = plt.subplots(figsize=(5, 5))
sns.scatterplot(data=_accepted, x="pilivar", y="dwell_time")

fig, ax = plt.subplots(figsize=(5, 5))
sns.scatterplot(data=_accepted, x="kb_sh", y="k_spawn")

# %%

def describe_std_ratio(a, b):
    print("{:.3f}/{:.3f} = {:.3f}".format(a,b,a/b))
describe_std_ratio(accparams["pdisp"].std(), params["pdisp"].std())
accparams[names + ["pdisp"]].corr()

# %%

accparams[names + _objectives].corr()

# %%
accparams
accparams['kmsd.mean'].mean()

# %%
params = sim4d["params"]
describe_std_ratio(accparams["ntaut.mean"].std(), params["ntaut.mean"].std())
accparams[names + ["ntaut.mean"]].corr()

# so are dwell_time/kspawn estimates anticorrelated then?
# -0.155 


# %%

# * params[names + ["qhat.estimate"]].corr()
# params[names + ["qhat.estimate"]].corr()
# params[names + ["qhat.estimate", "kmsd.mean", "ntaut.mean", "pdisp"]].corr()
# accparams[names + ["lvel.mean"]].corr()
# accparams[names + ["qhat.estimate", "kmsd.mean", "ntaut.mean", "pdisp"]].corr()

# sns.scatterplot(data=accparams, x="dwell_time", y="lvel.mean")
# sns.scatterplot(data=accparams, x="k_spawn", y="lvel.mean")

# %%
# * plotting the ntaut.mean / kmsd.mean anticorrelation
# fig, ax = plt.subplots(figsize=(5, 5))
# sns.scatterplot(data=accparams, x="ntaut.mean", y="kmsd.mean")

# %%
# * check simulation with similar qhat have similar kmsd.mean
# fix_qhat = params[params["qhat.estimate"].between(0.39,0.41)]
# scipy.stats.pearsonr(params["qhat.estimate"], params["kmsd.mean"])
# fig, ax = plt.subplots(figsize=(5,5))
# sns.scatterplot(data=params, x="qhat.estimate", y="kmsd.mean")
# sns.scatterplot(data=fix_qhat, x="qhat.estimate", y="kmsd.mean")
# todo pick a high kmsd, low qhat sample and a low kmsd high qhat sample to look at (?)


# %%
# plot scatter plots pdispj

fig, ax = plt.subplots(figsize=(5, 5))
sns.scatterplot(data=accparams, x="k_spawn", y="pdisp")

fig, ax = plt.subplots(figsize=(5, 5))
sns.scatterplot(data=accparams, x="kb_sh", y="pdisp")



# %%
# print 1d summary statistic distributions
_objectives = ["lvel.mean", "deviation.var", "qhat.estimate",
               "ahat.estimate", "nbound.mean", "ntaut.mean", "pdisp"]
abcimplement.plot_accepted_stats(
    sim4d["params"], _accepted, _objectives, simref)

# %%


scipy.stats.pearsonr(accparams['k_spawn'], accparams['ntaut.mean']), scipy.stats.pearsonr(accparams['dwell_time'], accparams['ntaut.mean'])
# scipy.stats.pearsonr(accparams['dwell_time'], accparams['k_spawn'])
# plt.scatter(accparams['dwell_time'], accparams['k_spawn'])

scipy.stats.pearsonr(accparams['kb_sh'], accparams['pdisp'])

# %%
# params['ntaut.mean'].std(), accparams['ntaut.mean'].std()
# params['ntaut.mean'].std()/ accparams['ntaut.mean'].std()
accparams['pdisp'].std(), params['pdisp'].std(), accparams['pdisp'].std()/params['pdisp'].std()
# params['ntaut.mean'].std(), accparams['ntaut.mean'].std()

# %%
# production quality distribution plots for ntaut.mean, pdisp.mean


def plot_kspawn_ntaut(ax):
    sns.scatterplot(data=accparams, x="k_spawn", y="ntaut.mean", ax=ax)
    ax.set_ylabel(r"$\langle N_\mathrm{taut} \rangle$", fontsize=20)
    ax.set_xlabel(r"$k_{\mathrm{spawn}}$", fontsize=20, labelpad=0)

def plot_tdwell_ntaut(ax):
    sns.scatterplot(data=accparams, x="dwell_time", y="ntaut.mean", ax=ax)
    ax.set_ylabel(r"$\langle N_\mathrm{taut} \rangle$", fontsize=20)
    ax.set_xlabel(r"$\tau_{\mathrm{dwell}}$", fontsize=20, labelpad=0, backgroundcolor="white")


with mpl.rc_context({"font.size": 20,"text.usetex":True}):
    fig, ax = plt.subplots(figsize=(5,5))
    defcolor = plt.rcParams['axes.prop_cycle'].by_key()['color']
    hstyle = {'stat':'density', 'alpha':0.4, 'element':'step'}
    sns.histplot(params["ntaut.mean"], ax=ax, color=defcolor[0], **hstyle)
    arr = accparams["ntaut.mean"].to_numpy()
    sns.histplot(arr, bins=8, ax=ax, color=defcolor[1],  **hstyle)
    ax.set_xlabel(r"$\langle N_{\mathrm{taut}} \rangle$", fontsize=24)
    ax.set_ylabel("Density",fontsize=24)
    plt.locator_params(axis='y', nbins=6)
    ax.legend(["all", "accepted"], fontsize=16)

    # inset
    inset = ax.inset_axes([0.35, 0.32, 0.4, .4])
    plot_kspawn_ntaut(inset)
    # inset.set_xlim((0.1,8.0))
    x = accparams["k_spawn"]
    inset.set_xticks([x.min(), x.max()])
    inset.set_xticks([0, 8])
    inset.set_xticklabels(['0', '8'], fontsize=16)
    # inset.set_yticks([])

    inset = ax.inset_axes([0.77, 0.32, 0.4, .4])
    # inset.set_xlim((0,np.pi/2))
    # why does mpl hide labesls when argument is empty list []?
    x = accparams["dwell_time"]
    plot_tdwell_ntaut(inset)
    inset.set_ylabel("")
    # inset.set_xticks([0])
    inset.set_xticks([0, 3.0])
    inset.set_xticklabels(['0', '3.0'], fontsize=16)
    inset.set_xlim([0,  3.0])
    inset.set_yticks([0])

pub.save_figure("ntaut_distrib_with_inset", notename)

#  %%

def plot_kspawn_pdisp(ax):
    sns.scatterplot(data=accparams, x="k_spawn", y="pdisp", ax=ax)
    ax.set_ylabel(r"$\langle \Delta x \rangle$", fontsize=20)
    ax.set_xlabel(r"$k_{\mathrm{spawn}}$", fontsize=20, labelpad=-6)

def plot_alpha_pdisp(ax):
    sns.scatterplot(data=accparams, x="kb_sh", y="pdisp", ax=ax)
    ax.set_ylabel(r"$\langle \Delta x \rangle$", fontsize=20)
    ax.set_xlabel(vret, fontsize=20, labelpad=0, backgroundcolor="white")
    # plt.locator_params(axis='x', nbins=3)
    
with mpl.rc_context({"font.size": 20, "text.usetex":True}):
    fig, ax = plt.subplots(figsize=(5,5))
    defcolor = plt.rcParams['axes.prop_cycle'].by_key()['color']
    hstyle = {'stat':'density', 'alpha':0.4, 'element':'step'}
    sns.histplot(params["pdisp"], ax=ax, color=defcolor[0], **hstyle)

    arr = accparams["pdisp"].to_numpy()
    sns.histplot(arr, bins=8, ax=ax, color=defcolor[1],  **hstyle)
    ax.set_xlabel(r"$\langle \Delta x \rangle$", fontsize=24)
    ax.set_ylabel("Density",fontsize=24)
    ax.set_ylim((0,18))
    ax.legend(["all", "accepted"], fontsize=16)

    # inset
    inset = ax.inset_axes([0.35, 0.32, 0.4, .4])
    plot_kspawn_pdisp(inset)
    # inset.set_xlim((0.1,8.0))
    x = accparams["k_spawn"]
    inset.set_ylabel(r"$\langle \Delta x \rangle$", labelpad=-20)
    # inset.set_xticks([x.min(), x.max()])
    inset.set_xticks([0, 8])
    inset.set_xticklabels(['0', '8'], fontsize=16)
    inset.set_yticks([0.04, 0.16])
    inset.tick_params(axis='y', which='major', labelsize=16)

    # inset.set_yticklabels(fontsize=16)


    inset = ax.inset_axes([0.77, 0.32, 0.4, .4])
    # inset.set_xlim((0,np.pi/2))
    # why does mpl hide labesls when argument is empty list []?
    x = accparams["kb_sh"]
    plot_alpha_pdisp(inset)
    inset.set_ylabel("")
    # inset.set_xticks([0])
    inset.set_xticks([0.1, 0.4])
    inset.set_xticklabels(['0.1', '0.4'], fontsize=16)
    inset.set_xlim([0,  0.4])
    inset.set_yticks([0])

pub.save_figure("pdisp_distrib_with_inset", notename)

#  %%
# production quality distribution plots for ntaut.mean, pdisp.mean   plt.locator_params(axis='y', nbins=6)


with mpl.rc_context({"font.size": 20}):
    fig, ax = plt.subplots(figsize=(4,4))
    plot_kspawn_pdisp(ax)

with mpl.rc_context({"font.size": 20}):
    fig, ax = plt.subplots(figsize=(4,4))
    plot_alpha_pdisp(ax)

#  %%

fig, ax = plt.subplots(figsize=(5, 5))
with mpl.rc_context({"font.size": 20}):
    sns.scatterplot(data=accparams, x="k_spawn", y="ntaut.mean", ax=ax)
    ax.set_ylabel(r"$\langle N_\mathrm{taut} \rangle$", fontsize=20)
    ax.set_xlabel(r"$k_{\mathrm{spawn}}$", fontsize=20, labelpad=0)

fig, ax = plt.subplots(figsize=(5, 5))
with mpl.rc_context({"font.size": 20}):
    sns.scatterplot(data=accparams, x="dwell_time", y="ntaut.mean", ax=ax)
    ax.set_ylabel(r"$\langle N_\mathrm{taut} \rangle$", fontsize=20)
    ax.set_xlabel(r"$\tau_{\mathrm{dwell}}$", fontsize=20, labelpad=0)



# fig, ax = plt.subplots(figsize=(5, 5))
# sns.scatterplot(data=accparams, x="anchor_angle_smoothing_fraction", y="pdisp")


# %%

# for obs in  ["pdisp"]:
for obs in _objectives + ["pdisp"]:

    _regdf = statdf[['uid'] + sim4d["problem"]["names"] + [obs]]
    _acc = abcimplement.rejection_abc(_regdf, [obs], statref, N)
    accparams = params.iloc[_acc.index]

    fig, ax = plt.subplots(figsize=(5, 5))
    X = accparams['k_spawn']
    Y = accparams['dwell_time']
    ax.scatter(X, Y)
    ax.set_title(obs)
    ax.set_xlabel('k_spawn')
    ax.set_ylabel('dwell_time')
    print(scipy.stats.pearsonr(X, Y))
    print()


# %%
# ----------------------------------------------------------------
# * for more analysis see simabc.py

# %%
# * ----------------------------------------------------------------
# * ----------------------------------------------------------------
# * copy using dstep=0.03 ()

# %%
# localjson = "light.json"
localjson = "no_coarse.json"

simtarget = "/home/dan/usb_twitching/run/825bd8f/target/t0"
with command.chdir(simtarget):
    ltarget = stats.load(localjson)
    args = parameters.thisread()
_simobjective = ['lvel.mean', 'deviation.var',
                 'qhat.estimate', 'ahat.estimate']
simref = {name: twutils.make_get(name)(ltarget) for name in _simobjective}
_interest = ['dwell_time', 'k_spawn', 'pilivar',
             'anchor_angle_smoothing_fraction']

simpar = {par: args.pget(par) for par in _interest}
simpar['anchor_angle_smoothing_fraction'] *= np.pi/2
simpar, simref
simref

# %%

light4d = {}
light4d["simdir"] = "/home/dan/usb_twitching/run/825bd8f/cluster/mc4d"
light4d["objectives"] = ['lvel.mean',
                         'deviation.var', 'qhat.estimate', 'ahat.estimate']
light4d = abcimplement.load_problem_simulation(
    light4d, local_stats_path=localjson)
light4d["problem"]

# %%
# ABC config
N = 50
print('{}/{}'.format(N, light4d["M"]))

# %%
# one statistic at a time
_objectives = light4d["objectives"]
light4d["params"] = light4d["data"].paramsdf(_objectives)
abcimplement.transform_anchor_parameter_data(light4d)
statdf, statref = abcimplement.regularise_stats(
    light4d["params"], simref, _objectives)

# %%
# all four simple metrics
warnings.filterwarnings("ignore")

_objectives = ["sim.t0.ks_statistic",
               "deviation.var", "qhat.estimate", "ahat.estimate"]
_objectives = ["lvel.mean", "deviation.var", "qhat.estimate", "ahat.estimate"]
# _objectives = ["lvel.mean", "deviation.var", "nbound.mean", "ahat.estimate"]
# _objectives = ["lvel.mean", "deviation.var", "nbound.mean"]
_objectives = ["lvel.mean", "deviation.var", "qhat.estimate", "ahat.estimate"]

_regdf = statdf[['uid'] + light4d["problem"]["names"] + _objectives]
_accepted = abcimplement.rejection_abc(_regdf, _objectives, statref, N)

df = abcimplement.tabulate_inference(
    light4d["problem"], _accepted, _objectives, simpar=simpar)
m_par = {k: v for k, v in zip(df["parameter"], df["MLE"])}

prime_accepted = _accepted
pub.save_dfhtml(df, "sim_crawlingabc", notename)
df

# %%
if True:
    fig, axes = abcimplement.perfectplot4d(
        light4d["problem"], _accepted, simpar=simpar, mpar=m_par)
    # pub.save_figure('new_sim_crawling_abc', notename, fig)

# %%
names = sim4d["problem"]["names"]
_accepted[names].corr()

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(_accepted["pilivar"], _accepted["dwell_time"])

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(_accepted["anchor_angle_smoothing_fraction"], _accepted["k_spawn"])


# %%
params = sim4d["params"]
accparams = params.iloc[_accepted.index]
accparams
