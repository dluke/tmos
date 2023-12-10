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
# towards publication of sobol data 
# note that sobolnote.py is a bit bloated at this point

# %% 
import copy
import os, sys
join = os.path.join
import numpy as np
import pandas as pd
import stats
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import SALib.analyze as sal
#
import rtw
import _fj
import txtdata
import fjanalysis
import twutils
import twanalyse
import plotutils
import sobol

import pili
from pili import publication
import pili.publication as pub
import seaborn as sns

# %% 
notedir = os.path.split(os.getcwd())
notename = "sobolpub"

# we will be loading several datasets, use dictionary for each dataset to keep global scope clean

# first dataset is 4d sobol that we use to check the search limits, particularly for dwell_time 
# but also we need to expand k_spawn to the right and \alpha to the left

# %% 
import collections

# combine loading the various data and meta data for a MC/sobol simulation with some 'error' checking
# todo: move to sobol module
def load_problem_simulation(psd): # problem simulation data
    _lookup = sobol.read_lookup(psd["simdir"])
    _problem = sobol.read_problem(psd["simdir"])
    # check failed
    types= {"failed_condition": object}
    Y = sobol.collect(["failed", "failed_condition"], targetdir=psd["simdir"], tdct=types)
    Y["failed"][np.isnan(Y["failed"])] = False
    failed = Y["failed"].astype(bool)
    # assert( all(Y["failed_condition"][failed] == 'step_condition') )
    failed_count = collections.Counter(Y["failed_condition"][failed])
    print("failed:", failed_count)
    nan_idx = np.argwhere(failed).ravel()
    mask = sobol.filter_missing(_problem, _lookup, nan_idx, calc_second_order=psd["second_order"])
    N = len(_lookup[0])
    # filter missing data
    for i, m in reversed(list(enumerate(mask))):
        if not m:
            del _lookup[0][i]
    print("filtered out {}/{} samples".format(N - len(_lookup[0]), N))
    psd["lookup"] = _lookup
    psd["problem"] = _problem
    _lduid = sobol.collect_lduid(psd["simdir"])
    psd["data"] = sobol.SobolData(_problem, _lookup, _lduid, psd["second_order"])
    print("loaded data from ", psd["simdir"])
    return psd

lim4d = {}
lim4d["simdir"] = "/home/dan/usb_twitching/run/5bfc8b9/cluster/sobol_tdwell"
lim4d["second_order"] = True
objectives = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate']
load_problem_simulation(lim4d)

print(lim4d["problem"])

# %%
Y = {name : lim4d["data"].get(name) for name in objectives}
S = lim4d["data"].compute(Y)
dftable1, dftableT = lim4d["data"].format(S)
dftableT 

# %%  [markdown]
# Sensitivity indices are similar but with key difference that 
# index(dwell_time, ahat.estimate) is now ~0.1, even larger than
# index(pilivar, ahat.estimate), suggesting that we should investigate
# very low dwell_time further to conclude whether it can explain the slow crawling behaviour
# and also that we should probably argue that dwell_time << 0.5 is unlikely
# and that will preserve the special charactistic of the ahat metric 
# and our plan for sequential ABC on crawling data.

# %%  
# On the the main purpose of the this note
# we will take some of the following parameters 
# [k_ret_on, k_ret_off, k_ext_on, args.pili.L, f_stall, ks, "k-retraction"]
# varying f_stall and ks is almost certainly redundant
# we will handle k_ret_on separately since this is the surface sensing parameter
# we should also refer to Koch et al. for the uncertainty in their parameters
part_problem = {
    "names" : ["k_ret_on", "k_ret_off", "k_ext_on", "Lp", "ks"],
    "bounds" : [
        [0.25,2.5],  # [0.4s, 4s]
        [0.05,0.20],  # [5.0, 20.0]
        [0.20,0.50],  # [2.0, 5.0]
        [1.0, 10.0],  # [Lp]
        [1.*10**3, 2*10**4],  # [ks]
    ]
}
# now setup these jobs and submit them one by one with, trail them against the first 5 parameters
# (N = 1024, calc_second_order = False)

# %%
# testdir  = "/home/dan/usb_twitching/run/825bd8f/cluster/sobol_plus"
testdir = "/home/dan/usb_twitching/run/latest/sobol"
_lookup = sobol.read_lookup(testdir)


# ! diff the lookup with the actual folders here

# len(_lookup[0]), 13 * 2048
len(_lookup[0]), 13 * 2048 # (11 + 2)

# %%
# !sobol
sten = {}

# sten["simdir"] = "/home/dan/usb_twitching/run/825bd8f/cluster/sobol"
# sten["simdir"] = "/home/dan/usb_twitching/run/825bd8f/cluster/sobol_plus"
# !latest
sten["simdir"] = "/home/dan/usb_twitching/run/latest/sobol"

sten["second_order"] = False
objectives = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate']

load_problem_simulation(sten)
print(sten["problem"])

# %%
n = len(sten["data"].lookup[0])

# %%
Y = {name : sten["data"].get(name) for name in objectives}
Y['lvel.mean'] = np.nan_to_num(Y['lvel.mean']) #! only 1 NAN value
S = sten["data"].compute(Y)
dftable1, dftableT = sten["data"].format(S)

# %%
# add column to the total
rowmax = dftableT[objectives].max(axis=1)
dftableT["max"] = rowmax

# %%

cm = sns.light_palette("green", as_cmap=True)
cmlist = [cm(i) for i in np.linspace(0, 1, 256, True)]
n = 6 
cmlist = cmlist[:80]
# cmlist[:n] = [cmlist[0] for _ in range(n)]
cmlist[:n] = [(1.0,1.0,1.0,1.0) for _ in range(n)]
cm = mpl.colors.LinearSegmentedColormap.from_list('modgreen', cmlist)
cm


# %%

sty = dftableT.style.hide_index()
sty.to_latex()
sty.background_gradient(cmap=cm)
# sty.text_gradient(cmap=cm)
publication.save_dflatex(sty, "allsobol", "sobolpub")
sty

# %%
# ( "{{{%s}}}" % "c"+str(i)+str(j) ).format(**{"c"+str(i)+str(j) : 'ffffff'}) = "{ffffff}"

# include min/max parameters

cfl = ["s", ".5g", ".5g", ".4f", ".4f", ".4f", ".4f", ".4f"]
body = pub.construct_template(cfl, 11)

# publication.write_tex(template, "all_sobol_part", notename)
head = r"{h0} & {h1} & {h2} & {h3} & {h4} & {h5} & {h6} & {h7} \\"

# load allsobol.template.template
tt = publication.load_template("allsobol.template")
template = tt.format(head=head, body=body)
# save to allsobol.template
publication.save_template("allsobol", template)

# %%

header = ["parameter", 'min', 'max', 
    r'$\mean{u}$', r'$Var({\theta_d})$', r'$ \hat{q} $', r'$ \hat{a} $', r"$\max{S_i^T}$"]

lpar = ['\\keoff', '\\tdwell', '\\kappa', '\\vret', '\\kspawn', '\\keon', 
    '\\kresample', '\\kroff', ' $L_p$', ' $E$', '\\fstall']

df = dftableT[objectives + ["max"]]
_hdct, _pdct, _cdct, _tdct, _vdct = publication.construct_dicts(df, 
        sten["problem"], objectives, cm, header=header, lpar=lpar)

# %%
template = publication.load_template("allsobol")



# %%
template = publication.load_template("allsobol")
# need to cut first and last lines
_lines = template.strip().split('\n')
template = r'\n'.join(_lines[1:-1])
filled = template.format(**_hdct, **_pdct, **_cdct, **_tdct, **_vdct)
filled = _lines[0]+r'\n'+filled+r'\n'+_lines[-1]
publication.write_tex(filled, "allsobol_gen", notename)
# filled

# %%
