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



INDIR = "Profiles_pop3/Normal_profiles"
OUTDIR = "Profiles_pop3/Normal_graphs_lin"
OUTFILE = "sedov_taylor_distances"
FIELDS = ['density', 'temperature', 'velocity_magnitude', 'metallicity3']
RHO_BKGRD = unyt_quantity(1e-25, 'g/cm**3')


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = glob.glob(os.path.join(INDIR, "DD*_profile_weight_field_cell_mass.h5"))
datasets.sort(key = lambda x: get_time_z(os.path.basename(x), star_type)[0])
data_distance = [unyt_quantity(0.0, 'pc') for i in range (len(datasets))]
e_factor = (get_sn_energy(star_type) / RHO_BKGRD)**(1/5)
time_offset = get_lifetime_offset(star_type)
theory_distance = [unyt_quantity(0.0, 'pc') for i in range (len(datasets))]
times = [unyt_quantity(0.0, 'Myr') for i in range (len(datasets))]

for field in FIELDS:
    for i, dataset in enumerate(datasets):
        
        dump_name = os.path.basename(dataset)
        ds = yt.load(dataset)
        field_dict = get_field_dict(field)
        used = ds.profile.used
        radii = ds.profile.x[used].to('pc')

        if (field == 'temperature'):
            temps = ds.data[('data', 'temperature')][used]
            grads = [((temps[i+1] - temps[i]) / temps[i]) for i in range (len(temps) - 1)]
            index_max = np.argmin(transpose_unyt(grads)) + 1
            data_distance[i] = radii[index_max].to('pc')
            times[i] = get_time_z(dump_name, star_type)[0]
            theory_distance[i] = e_factor * (times[i] - time_offset)**(2/5)

            
        plt.figure()
        plt.title(get_title(dump_name, star_type))
        plt.plot(radii, ds.data[('data', field)][used])
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
        plt.xlabel("Radius (pc)")
        #plt.xlim(1e-1, 1e4)
        plt.xlim(0, 400)
        #plt.xscale('log')
        if (field_dict['log'] == True):
            plt.yscale('log')
        if (field == 'temperature'):
            plt.axvline(x= radii[index_max])
        plt.ylim(field_dict['limits'])
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_{field}.png"))
        plt.close()


'''
data_distance = transpose_unyt(data_distance).to('pc')
theory_distance = transpose_unyt(theory_distance).to('pc')
times = transpose_unyt(times)
data_distance.write_hdf5(f"{OUTFILE}.h5", group_name="data")
theory_distance.write_hdf5(f"{OUTFILE}.h5", group_name="theory")
times.write_hdf5(f"{OUTFILE}.h5", group_name="times")
'''
