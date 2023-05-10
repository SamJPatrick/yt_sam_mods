import os
import yt
import sys
import glob
import numpy as np
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc


FILENAME = "sedov_taylor_distances.h5"
OUTFILE = "sedov_taylor_distances.png"


data_distance = transpose_unyt(data_distance).to('pc')
theory_distance = transpose_unyt(theory_distance).to('pc')
times = transpose_unyt(times)

df = h5py.File(FILENAME, 'r')
times = list(df['times']['array_data'])
data_distance = list(df['data']['array_data'])
theory_distance = list(df['theory']['array_data'])
plt.figure()
plt.title("Blast wave radius")
plt.plot(times, data_distance, label= 'data', color='orange')
plt.plot(times, theory_distance, label= 'theory', color='blue')
plt.ylabel("Radius (pc)")
plt.xlabel("Time (Myr)")
plt.xlim(0.0, 30.0)
plt.ylim(0.0, 1000.0)
plt.legend()
plt.savefig(f"{OUTFILE}.png")
plt.close()
