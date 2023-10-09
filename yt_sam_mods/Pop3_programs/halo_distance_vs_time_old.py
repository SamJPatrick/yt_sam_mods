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

RADIUS_CUTOFF = unyt_quantity(250.0, 'pc')
DIVISOR = 3

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
    x_front = [0.0 for j in range (len(times_ray))]
    for j in range (len(times_ray)):
        df = h5py.File(ray_files[-j-1], 'r')
        distances_ray = unyt_array(df["distances/array_data"], 'pc')
        index_cutoff = np.argwhere((distances_ray - RADIUS_CUTOFF) > 0)[0].item()
        temp_trans = [unyt_array(np.round(x, decimals=3), 'K') for x in zip(*[df[f'temperature_{n}/array_data'] for n in range (NUM_RAYS)])]
        means_raw = [unyt_quantity(np.round(np.mean(arr), decimals=3), 'K') for arr in temp_trans]
        ranges = transpose_unyt([np.ptp(arr) for arr in temp_trans])
        #means = [np.mean([x for x in arr if (x <= (means_raw[k] + ranges[k]/2) and x >= (means_raw[k] - ranges[k]/2))]) for k, arr in enumerate(temp_trans)]
        means = [0.0 for k in range (len(temp_trans))]
        for k, arr in enumerate(temp_trans):
            mean_arr = []
            for x in arr:
                if ((x <= means_raw[k] + ranges[k]/DIVISOR) and (x >= means_raw[k] - ranges[k]/DIVISOR)):
                    mean_arr.append(x)
            if (len(mean_arr) == 0):
                print("DANGER, DANGER WILL WILSON!!!!!!!!")
                print(df.file)
                print(k)
                print(distances_ray[k])
                print(means_raw[k])
                print(ranges[k])
                print(arr)
                means[k] = means_raw[k]
            else :
                means[k] = np.mean(transpose_unyt(mean_arr))
        grads = [((means[k+1] - means[k]) / means[k]) for k in range (len(means) - 1)]
        indices = np.argsort(grads[:index_cutoff])
        for index in indices:
            if (distances_ray[index] > x_font[-j]):
                blah;
        x_front[-j-1] = df[f"distances/array_data"][indices[0]+1]
    plt.plot(times_ray, x_front, color= COLORS[i], linestyle='--')

plt.title("Halo & shock distances vs time")
plt.legend(loc= 'upper right')
plt.ylabel("Distance (pc)")
plt.xlabel("Time (Myr)")
plt.ylim(0, 600)
plt.xlim(-5, 30)
plt.axvline(x=0.0, color='black')
plt.savefig(os.path.join(OUTDIR, "distance_vs_time_test.png"))
plt.close()
