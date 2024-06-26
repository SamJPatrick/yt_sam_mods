import os
import sys
import csv
import glob
import h5py
from unyt import unyt_array
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt


DISTANCE_FILE = "ray_distances.txt"
SIM_FILE = "simulation.h5"

RAY_DIR = "Sobolev/Ray_profiles"
OUTDIR = '.'
STAR_DIRS = ['CCSN', 'HN', 'PISN']
COLORS = ['green', 'blue', 'red']

NUM_RAYS = 20
Z0 = unyt_quantity(1e-7, 'Zsun')



plt.figure()
halo_files = [os.path.join(star_dir, DISTANCE_FILE) for star_dir in STAR_DIRS]
sim_files = [os.path.join(star_dir, SIM_FILE) for star_dir in STAR_DIRS]
for i, halo_file in enumerate(halo_files):
    
    star_type = STAR_DIRS[i].lower()
    with open(halo_file, newline='\n') as myfile:
        reader = list(csv.reader(myfile, delimiter= '\t'))
    dumps, distances = zip(*[entry for entry in reader])
    distances_halo = unyt_array([float(distance) for distance in distances], 'pc')
    times_halo = transpose_unyt([get_time_z(filename, star_type, sim_files[i])[0] for filename in dumps])
    plt.plot(times_halo, distances_halo, color = COLORS[i], label= star_type)

    ray_files = sorted(glob.glob(os.path.join(STAR_DIRS[i], RAY_DIR, "DD*_packed.h5")))
    times_ray = transpose_unyt([get_time_z(os.path.basename(filename), star_type, sim_files[i])[0] for filename in ray_files])
    times_ray = times_ray[times_ray > get_lifetime_offset(star_type)]
    ray_files = ray_files[len(ray_files)-len(times_ray):]
    
    x_front = [0.0 for j in range (len(times_ray))]
    for j in range (len(times_ray)):
        df = h5py.File(ray_files[j], 'r')
        distances_ray = unyt_array(df["distances/array_data"], 'pc')
        mean_metals = unyt_array([np.mean(x) for x in zip(*[df[f'metallicity3_{n}/array_data'] for n in range (NUM_RAYS)])], 'Zsun')
        indices = np.where(mean_metals < Z0)[0]
        if (len(indices) == 0):
            x_front[j] = x_front[j-1]
            continue
        for index in indices:
            if (j==0 or distances_ray[index] >= x_front[j-1]):
                x_front[j] = distances_ray[index]
                break
        if (x_front[j] == 0.0):
            x_front[j] = x_front[j-1]
    plt.plot(times_ray, x_front, color= COLORS[i], linestyle='--')
    print(f"Blast wave distances for {STAR_DIRS[i]} is {x_front}")

plt.title("Halo & shock distances vs time")
plt.legend(loc= 'upper left')
plt.ylabel("Distance (pc)")
plt.xlabel("Time (Myr)")
plt.ylim(0, 600)
plt.xlim(-5, 60)
plt.axvline(x=0.0, color='black')
plt.savefig(os.path.join(OUTDIR, "distance_vs_time_test.png"))
plt.close()
