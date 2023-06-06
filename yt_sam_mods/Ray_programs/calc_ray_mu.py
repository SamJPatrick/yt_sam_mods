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
OUTDIR = "Sobolev/Ray_graphs"
#STD_ANGLE = 0.013784005317407418
WIN = 1



try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

ray_files = sorted(glob.glob(os.path.join(INDIR, "DD*_packed.h5")))
for ray_file in ray_files:
    
    df = h5py.File(ray_file, 'r')
    distances = unyt_array(df['distances'][:], 'pc')
    arr_len = len(distances)
    std_psi = transpose_unyt([np.std(transpose_unyt(vels)) for vels in \
                              zip(*[unyt_array(df[f'velocity_{n}'][:], 'km/s') for n in range (NUM_RAYS)])]) / np.std(df['angles'][:])
    std_rad = transpose_unyt([np.mean(transpose_unyt([np.std(unyt_array(df[f'velocity_{n}'][:][i-WIN:i+WIN+1], 'km/s')) / (distances[i-WIN] - distances[i+WIN]) \
                                       for n in range (NUM_RAYS)])) for i in range (WIN, arr_len - WIN)])
    velocities = transpose_unyt([np.mean(transpose_unyt(vels)) for vels in zip(*[unyt_array(df[f'velocity_{n}'][:], 'km/s') for n in range (NUM_RAYS)])])
    pressures = transpose_unyt([np.mean(transpose_unyt(press)) for press in zip(*[unyt_array(df[f'pressure_{n}'][:], 'Pa') for n in range (NUM_RAYS)])])
    velocities_diff = transpose_unyt([(velocities[i+1] - velocities[i]) / (distances[i+1] - distances[i]) for i in range (WIN, arr_len - WIN)])
    mu = (0.5 * pressures[WIN:arr_len-WIN] * np.abs(velocities_diff) / [(6 / distances[WIN:arr_len-WIN]**2) * std_psi[WIN:arr_len-WIN]**2 + std_rad**2])[0]

    plt.figure()
    dump_name = os.path.basename(ray_file)
    plt.title(get_title(dump_name, star_type))
    plt.plot(distances[WIN:arr_len-WIN], mu)
    plt.ylabel(f"$\mu (Pa s)$")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.yscale('log')
    plt.ylim(1e-13, 1e-2)
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_mu.png"))
    plt.close()
