# %% [markdown]
# a reviewer doesn't like that we selected a trajectories by hand 
# so lets fix that for him

# %% 
import os
import sys
import json
join = lambda *x: os.path.abspath(os.path.join(*x))
import numpy as np
norm = np.linalg.norm
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import shapeplot
import filesystem

import _fj
import plotutils
import command
import twanalyse
import fjanalysis
import twanimation
import pili
import stats

import pandas as pd
import seaborn as sns
import pili.publication as pub

import pili.support as support
from pili.support import make_get, get_array

from tqdm import tqdm

# %% 
df = pd.read_pickle("classification.pkl")
subset = _fj.load_subsets()

# trajectries showing circle anomaly are clearly incorrectly tracked or not twitching
# so lets just remove them from consideration
valid = ~np.logical_or.reduce([np.isnan(df[key]) for key in df]) & ~df['circle_anomaly']
len(df)-valid.sum() # removed 161

# %% 
np.nanquantile(df['lvel'], 0.1), np.nanquantile(df['lvel'], 0.9)
# np.nanquantile(df['lvel'], 0.25), np.nanquantile(df['lvel'], 0.75)

# %% 
top = subset["top"]
walking = subset["walking"]

# our goal here is to replace the 'top' and 'walking' subsets
# 'top' comes from the top 100 mean velocity crawling tracks
# 'walking' need to filter the static trajectories out of the walking set 

require_duration = 100
require_nsteps = 100
require_displacement = 3 
# require_walking_displacement = 1
require_aspect = 1.6

short = np.logical_or.reduce([~valid, df["duration"] < require_duration, df["nsteps"] < require_nsteps])
high_aspect = df["min_aspect"] > require_aspect
high_displacement = df["displacement"] > require_displacement

crawling = np.logical_and.reduce([~short, high_aspect, high_displacement])
walking = np.logical_and.reduce([~short, ~high_aspect, high_displacement])
trapped = np.logical_and.reduce([~short, high_aspect, ~high_displacement])
wtrapped = np.logical_and.reduce([~short, ~high_aspect, ~high_displacement])
horizontal = np.logical_or(trapped, crawling)

print("identified {} short trajectories".format(np.sum(short)))
print("identified {} crawling trajectories".format(np.sum(crawling)))
print("identified {} walking trajectories".format(np.sum(walking)))

# low_aspect_idx = df["min_aspect"].to_numpy() < 1.6
# high_aspect_idx = low_aspect_idx[~


# %% 
# write these trajectory lists to new 'new_top' and 'new_walking' subsets
tdf = df.loc[crawling].sort_values("lvel").iloc[-100:]
wdf = df.loc[walking]
len(wdf)


# %% 

cache_dir = "plos_classification/"
filesystem.safemkdir(cache_dir)
np.savetxt(join(cache_dir, "new_crawling.npy"), tdf.index.to_numpy())
np.savetxt(join(cache_dir, "new_walking.npy"), wdf.index.to_numpy())


# %%
