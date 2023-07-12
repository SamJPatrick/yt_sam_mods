import os
import sys
import csv
from unyt import unyt_array
from yt.extensions.sam_mods.graph_funcs import get_time_z
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt


DISTANCE_FILE = "ray_distances.txt"
OUTDIR = '.'


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


with open(DISTANCE_FILE, newline='\n') as myfile:
    reader = list(csv.reader(myfile, delimiter= '\t'))
dumps, distances = zip(*[entry for entry in reader])
distances = unyt_array([float(distance) for distance in distances], 'pc')
times = transpose_unyt([get_time_z(filename, star_type)[0] for filename in dumps])

plt.figure()
plt.title("Halo distance vs time")
plt.plot(times, distances)
plt.ylabel("Inter-halo distance (pc)")
plt.xlabel("Time (Myr)")
plt.ylim(0, 1200)
plt.xlim(0, 140)
plt.savefig(os.path.join(OUTDIR, "halo_distance_vs_time.png"))
plt.close()
