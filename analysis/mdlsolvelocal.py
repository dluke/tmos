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
# similar to mdlpartition.py
# try to code some adjustments to the local solver

# %% 
import os
import json
import pickle
import numpy as np
join = lambda *x: os.path.abspath(os.path.join(*x))
norm = np.linalg.norm
import pandas as pd
from glob import glob

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import sctml
import sctml.publication as pub
print("writing figures to", pub.writedir)

import pili
import mdl
import pwlpartition

# %% 
target = "/home/dan/usb_twitching/sparseml/run/synthetic/inter_vary_sigma/sigma_00.1000/solver"
solver = pwlpartition.Solver.load_state(target)

fig, ax = plt.subplots(figsize=(20,20))
pwlpartition.simple_model_plot(ax, solver.partition.model, solver.partition.data)

# %% 
# create at 12
solver.partition.create_at(12)

fig, ax = plt.subplots(figsize=(20,20))
pwlpartition.simple_model_plot(ax, solver.partition.model, solver.partition.data)

# %% 
solver.percolate_at(12)

fig, ax = plt.subplots(figsize=(20,20))
ax.set_xlim(5,15)
pwlpartition.simple_model_plot(ax, solver.partition.model, solver.partition.data)


# %% 
solver.binary_percolate_at(12)

fig, ax = plt.subplots(figsize=(20,20))
ax.set_xlim(5,15)
pwlpartition.simple_model_plot(ax, solver.partition.model, solver.partition.data)

# %% 
for index in range(solver.partition.model.M-1):
    print('binary percolate at index',  index)
    solver.binary_percolate_at(index)
    
# %% 
fig, ax = plt.subplots(figsize=(40,40))
# ax.set_xlim(5,15)
pwlpartition.simple_model_plot(ax, solver.partition.model, solver.partition.data)

