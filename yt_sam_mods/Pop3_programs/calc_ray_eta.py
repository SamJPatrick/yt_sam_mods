import glob
import h5py
import numpy as np
import sys
import os

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.graph_funcs import *




NUM_RAYS = 10
INDIR = "Sobolev/Ray_profiles/Packed_rays"
OUTDIR = "Sobolev"


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

ray_files = sorted(glob.glob(os.path.join(INDIR, "DD*_packed.h5")))
eta = unyt_array([0.0 for i in range (len(ray_files))], 'Pa')
for i, ray_file in enumerate(ray_files):
    
    df = h5py.File(ray_file, 'r')
    distances = unyt_array(df['distances'][:], 'pc')
    arr_len = len(distances)
    densities = transpose_unyt([np.mean(transpose_unyt(dens)) for dens in zip(*[unyt_array(df[f'density_{n}'][:], 'g/cm**3') for n in range (NUM_RAYS)])])
    velocities = transpose_unyt([np.mean(transpose_unyt(vels)) for vels in zip(*[unyt_array(df[f'velocity_{n}'][:], 'km/s') for n in range (NUM_RAYS)])])
    pressures = transpose_unyt([np.mean(transpose_unyt(press)) for press in zip(*[unyt_array(df[f'pressure_{n}'][:], 'Pa') for n in range (NUM_RAYS)])])
    velocities_diff = transpose_unyt([(velocities[i+1] - velocities[i]) / (distances[i+1] - distances[i]) for i in range (arr_len - 1)])

    index_peak = np.argmax(np.nan_to_num(densities))
    index_final = np.argwhere(np.nan_to_num(densities[index_peak:]) < 0.1 * np.max(np.nan_to_num(densities[index_peak:])))[0].item() + index_peak
    eta_ray = unyt_quantity(0.0, 'Pa')
    for j in range (index_peak, index_final):
        eta_ray += abs(pressures[j] * velocities_diff[j] * (distances[j+1] - distances[j]) / velocities[j])
    eta[i] = eta_ray


#datasets.sort(key = lambda x: get_time_z(os.path.basename(x), star_type)[0])
plt.figure()
plt.title("Eta vs time PISN")
plt.plot([get_time_z(os.path.basename(x), star_type)[0] for x in ray_files], eta)
plt.ylabel(f"$\eta (J)$")
plt.xlabel("Time (Myr)")
plt.xlim(0.0, 130.0)
plt.yscale('log')
#plt.ylim(1e-13, 1e-2)
plt.savefig(os.path.join(OUTDIR, "eta_vs_time.png"))
plt.close()
