import glob
import h5py
import numpy as np
import sys
import os
import csv

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from unyt import unyt_quantity, unyt_array
from unyt import kb, mh
import yt
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.graph_funcs import *




NUM_RAYS = 10
WIN = 1

INDIR_RAYS = "Sobolev/Ray_profiles"
INDIR_PROFS = "Profiles_pop3/Normal_profiles"
DISTANCE_FILE = "ray_distances.txt"
OUTDIR_MU = "Sobolev/Viscous_params/Mu"
OUTDIR_NU = "Sobolev/Viscous_params/Nu"
OUTDIR_PHI = "Sobolev/Viscous_params/Phi"

MU = 1.6
CONV_CONST = kb / (MU * mh)


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


with open(DISTANCE_FILE, newline='\n') as myfile:
    reader = list(csv.reader(myfile, delimiter= '\t'))
dumps, distances = zip(*[entry for entry in reader])
halo_distances = unyt_array([float(distance) for distance in distances], 'pc')

prof_files = sorted(glob.glob(os.path.join(INDIR_PROFS, "DD*_profile_weight_field_cell_mass.h5")))
ray_files = sorted(glob.glob(os.path.join(INDIR_RAYS, "DD*_packed.h5")))
assert len(prof_files) == len(ray_files), "Error, unequal number of datasets found"
for i, ray_file in enumerate(ray_files):
    
    df = h5py.File(ray_file, 'r')
    dump_name = os.path.basename(ray_file)
    distances = unyt_array(df['distances/array_data'], 'pc')
    arr_len = len(distances)
    
    if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
        grad_field = "temperature"
        line_color = "red"
        df_prof = yt.load(prof_files[i])
        used = df_prof.profile.used
        radii = df_prof.profile.x[used].to('pc')
        temps = df_prof.data[('data', 'temperature')][used].to('K')
        grads = [((temps[i+1] - temps[i]) / temps[i]) for i in range (len(temps) - 1)]
        radius_max_prof = radii[np.argmin(grads) + 1].to("pc")
        yt.mylog.info(f"Profile radius of ST in {dump_name} is {radius_max_prof}")
    else :
        grad_field = "El_fraction"
        line_color = "green"
    means = [np.mean(x) for x in zip(*[df[f'{grad_field}_{n}/array_data'] for n in range (NUM_RAYS)])]
    grads = [((means[i+1] - means[i]) / means[i]) for i in range (len(means) - 1)]
    index_max = np.argmin(grads) + 1

    v_para = [unyt_array(df[f'velocity_para_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)]
    v_norm = [unyt_array(df[f'velocity_norm_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)]
    v_para_trans = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_para_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    v_norm_trans = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_norm_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    del_phi = np.std(df['angle_list/array_data'])

    a = transpose_unyt([unyt_quantity(np.std([np.mean(vel[i-WIN:i+WIN+1]) for vel in v_norm]), 'km/s') for i in range (WIN, arr_len-WIN)]) / del_phi
    b = transpose_unyt([(unyt_quantity(np.std([np.mean(vels) for vels in v_para_trans[i]]), 'km/s') / (distances[i+WIN] - distances[i-WIN])) for i in range (WIN, arr_len-WIN)])
    c = transpose_unyt([(unyt_quantity(np.std([np.mean(vels) for vels in v_norm_trans[i]]), 'km/s') / (distances[i+WIN] - distances[i-WIN])) for i in range (WIN, arr_len-WIN)]) / del_phi
    d = transpose_unyt([(unyt_quantity(np.std([np.mean(vel[i-WIN:i+WIN+1]) for vel in v_para]), 'km/s') / (distances[i+WIN] - distances[i-WIN])) for i in range (WIN, arr_len-WIN)]) / del_phi
    e = transpose_unyt([np.mean(vels) for vels in v_norm_trans[WIN:-WIN]])
    inv2_radii = transpose_unyt([1 / radius**2 for radius in distances[WIN:-WIN]])
    denom = 6 * inv2_radii * a**2 + b**2 + (c + d)**2 + 1.5 * inv2_radii * e**2

    velocities = transpose_unyt([np.mean(vels) for vels in v_para_trans])
    velocities_diff = transpose_unyt([(velocities[i+1] - velocities[i]) / (distances[i+1] - distances[i]) for i in range (WIN, arr_len-WIN)])
    temperatures = transpose_unyt([np.mean(transpose_unyt(temps)) for temps in zip(*[unyt_array(df[f'temperature_{n}/array_data'], 'K') for n in range (NUM_RAYS)])])
    temperatures_diff = transpose_unyt([(temperatures[i+1] - temperatures[i]) / (distances[i+1] - distances[i]) for i in range (WIN, arr_len-WIN)])
    cooling_times = transpose_unyt([np.mean(transpose_unyt(times)) for times in zip(*[unyt_array(df[f'cooling_time_{n}/array_data'], 's') for n in range (NUM_RAYS)])])
    densities = transpose_unyt([np.mean(transpose_unyt(dens)) for dens in zip(*[unyt_array(df[f'density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])
    numer = (temperatures[WIN:-WIN] / cooling_times[WIN:-WIN]) + (velocities[WIN:-WIN] * temperatures_diff) - (temperatures[WIN:-WIN] * velocities_diff)

    shear = denom.to('s**(-2)')
    nu = (0.5 * CONV_CONST * (numer / denom)).to('m**2/s')
    mu = (densities[WIN:-WIN] * nu).to('Pa*s')

    plt.figure()
    plt.title(get_title(dump_name, star_type))
    plt.plot(distances[WIN:-WIN], shear)
    plt.ylabel(r"$\Phi (s^{-2})$")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.yscale('log')
    plt.ylim(1e-17, 1e-10)
    plt.axvline(x= distances[index_max], color= line_color)
    plt.axvline(x= halo_distances[i], color='blue')
    if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
        plt.axvline(x= radius_max_prof, color='red', linestyle=':')
    plt.savefig(os.path.join(OUTDIR_PHI, f"DD{get_dump_num(dump_name)}_phi.png"))
    plt.close()

    plt.figure()
    plt.title(get_title(dump_name, star_type))
    plt.plot(distances[WIN:-WIN], mu)
    plt.ylabel(r"$\mu (Pa ~s)$")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.yscale('log')
    plt.ylim(1e-10, 1e4)
    plt.axvline(x= distances[index_max], color= line_color)
    plt.axvline(x= halo_distances[i], color='blue')
    if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
        plt.axvline(x= radius_max_prof, color='red', linestyle=':')
    plt.savefig(os.path.join(OUTDIR_MU, f"DD{get_dump_num(dump_name)}_mu.png"))
    plt.close()

    plt.figure()
    plt.title(get_title(dump_name, star_type))
    plt.plot(distances[WIN:-WIN], nu)
    plt.ylabel(r"$\nu (m^{2} ~s^{-1})$")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.yscale('log')
    plt.ylim(1e15, 1e25)
    plt.axvline(x= distances[index_max], color= line_color)
    plt.axvline(x= halo_distances[i], color='blue')
    if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
        plt.axvline(x= radius_max_prof, color='red', linestyle=':')
    plt.savefig(os.path.join(OUTDIR_NU, f"DD{get_dump_num(dump_name)}_nu.png"))
    plt.close()
