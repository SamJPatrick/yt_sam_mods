import os
import h5py
import glob
import yt
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_dump_num, get_time_z
from yt.extensions.sam_mods.misc import transpose_unyt



NORM_FOLDER = "Profiles/Normal_profiles"


filenames = glob.glob(os.path.join(os.getcwd(), NORM_FOLDER, "DD*_profile_weight_field_cell_mass.h5"))
filenames.sort(key= lambda x: int(get_dump_num(os.path.basename(x))))
mdot = [0.0] * len(filenames)
times = [0.0] * len(filenames)
for i, fn in enumerate(filenames):
    df = yt.load(fn)
    used = df.profile.used
    times[i] = get_time_z(os.path.basename(fn))[0]
    mdot[i] = df.profile[('data', 'accretion_rate')][used][-1]

times = transpose_unyt(times)
mdot =  transpose_unyt(mdot)
plt.figure()
plt.title("Accretion rate at virial radius vs time")
plt.plot(times, mdot)
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\dot{M} ~(M_{\odot} ~yr^{-1})$")
plt.axhline(y=0.0, color='black')
plt.ylim(-0.05, 0.25)
plt.savefig("virial_accretion_vs_time.png")
plt.close()
