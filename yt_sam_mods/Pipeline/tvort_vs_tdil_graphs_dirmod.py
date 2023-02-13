#########################################################################
# MODIFIED VERSION OF TVORT_VS_TDIR_GRAPHS.PY WHICH LOOKS UP WHETHER IMAGE IS ALRREADY IN TARGET DIRECTORY #
# POSSIBLY TO BE IMPLIMENTED IN FUTURE VERSIONS OF OTHER SCRIPTS #
##########################################################################

import os
import glob
import yt
import sys
from unyt import unyt_array

import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.misc import transpose_unyt



NORMDIR = "Profiles/Normal_profiles"
VELDIR = "Profiles/Velocity_profiles"
OUTDIR = "Profiles/Tvort_vs_tdil_graphs"


try:
    star_mode = sys.argv[1]
except IndexError:
    star_mode = ""
    pass

dumps_pre = glob.glob(os.path.join(os.getcwd(), OUTDIR, "DD*_tdil_vs_tvot.png"))


fn_norm = glob.glob(os.path.join(os.getcwd(), NORMDIR, "DD*_profile_weight_field_cell_mass.h5"))
fn_norm.sort(key= lambda x: int(get_dump_num(os.path.basename(x))))
fn_vel = glob.glob(os.path.join(os.getcwd(), VELDIR, "DD*_1D_profile_radius_cell_mass.h5"))
fn_vel.sort(key= lambda x: int(get_dump_num(os.path.basename(x))))
assert len(fn_norm) == len(fn_vel), "Error, directories must have equal number of nodes"

for i in range(len(fn_norm)):

    name = os.path.basename(fn_norm[i])
    _dump = f"DD{get_dump_num(name)}" + "_tdil_vs_tvort.png"
    if (_dump in dumps_pre):
       continue
    title = get_title(name, star_mode)
    dump = f"DD{get_dump_num(name)}"
    normdf = yt.load(fn_norm[i])
    veldf = yt.load(fn_vel[i])
    
    used = normdf.profile.used
    radii = normdf.profile.x[used].to("pc")
    log_sig_z = np.log(normdf.profile.standard_deviation[('data', 'metallicity3')][used])
    vr = veldf.profile[('data', 'velocity_spherical_phi')][used]
    t_dil = [np.abs(1  / (vr[j] * ((log_sig_z[j+1] - log_sig_z[j]) / (radii[j+1] - radii[j])))) for j in range(len(radii)-1)]
    t_dil = transpose_unyt(t_dil).to('Myr')
    t_vort = veldf.profile[('data', 'vortical_time')][used].to('Myr')

    plt.figure()
    plt.title(title)
    plt.plot(radii[:-1], t_dil, color='blue', label='dilution')
    plt.plot(radii, t_vort, color='orange', label='vortical')
    plt.xlabel("Radius (pc)")
    plt.xscale('log')
    plt.xlim(1e-2, 5e2)
    plt.ylabel(r"$t ~(Myr)$")
    plt.yscale('log')
    plt.ylim(5e-3, 5e4)
    plt.legend(loc= 'upper left')
    plt.savefig(os.path.join(OUTDIR, f"{dump}_tdil_vs_tvort.png"))
    plt.close()
