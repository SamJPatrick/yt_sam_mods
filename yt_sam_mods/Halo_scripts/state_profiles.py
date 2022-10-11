import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import return_sphere, align_sphere, transpose_unyt, modify_grackle
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, node_profile, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc

from grid_figure import GridFigure




def create_state_profiles(node, indir_vel='.', infile_be='.', outdir='.'):
   
    dsfn = node.ds_filename.split('/')[-1]
    filename = os.path.join(indir_vel, f"{str(dsfn)}_1D_profile_radius_cell_mass.h5")
    df = yt.load(filename)
    used = df.profile.used
    x_data = df.profile.x[used].to('pc')
    time = df.current_time.to('Myr')

    be_file = h5py.File(infile_be, 'r')
    be_times = unyt_array(list(be_file['time/array_data']), 'Myr')
    index, t_diff = min([(i, time - be_times[i]) for i in range(len(be_times))], key=lambda x: x[1]) 
    be_mass = be_file['mass_norm/array_data'][index]
    be_radius = be_file['radius_norm/array_data'][index]

    profile_dict = {}
    fields = ["turbulent_sound_crossing_time", "sound_crossing_time", "total_dynamical_time", \
              "cooling_time", "vortical_time"]
    for field in fields:
        if (field == "sound_crossing_time"):
            cs = df.profile[('data', 'sound_speed')][used]
            y_sc = (2 * x_data / cs)
            profile_dict[field] = y_sc.to('yr')
        elif (field == "turbulent_sound_crossing_time"):
            #cs = df.profile[('data', 'sound_speed')][used]
            vt = df.profile.standard_deviation[('data', 'velocity_magnitude')][used]
            #v = np.sqrt(cs**2 + vt**2)
            y_sct = (2 * x_data / vt)
            profile_dict[field] = y_sct.to('yr')
        else:
            profile_dict[field] = df.profile[('data', field)][used].to('yr')

    cell_mass = df.profile.weight[used].in_units('Msun')
    state_dict = {'frag':unyt_array([0.0 for i in range (len(x_data))], 'Msun'),
                  'support_turb':unyt_array([0.0 for i in range (len(x_data))], 'Msun'),
                  'collapse':unyt_array([0.0 for i in range (len(x_data))], 'Msun'),
                  'support_pressure':unyt_array([0.0 for i in range (len(x_data))], 'Msun')}
    for i in range (len(x_data) - 1):
        if (profile_dict['cooling_time'][i] < profile_dict['total_dynamical_time'][i]):
            state_dict['frag'][i] = cell_mass[i]
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['turbulent_sound_crossing_time'][i] and \
              profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            state_dict['collapse'][i] = cell_mass[i]
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            state_dict['support_turb'][i] = cell_mass[i]
        else :
            state_dict['support_pressure'][i] = cell_mass[i]

    gas_mass = unyt_array(np.zeros(len(state_dict.keys())), 'Msun')
    gas_frac = unyt_array(np.zeros(len(state_dict.keys())), '')
    mean_radius = unyt_array(np.zeros(len(state_dict.keys())), 'pc')
    radius_frac = unyt_array(np.zeros(len(state_dict.keys())), '')
    for i, state in enumerate(state_dict.keys()):
        gas_mass[i] = sum(state_dict[state])
        #gas_frac[i] = gas_mass[i] / sum(cell_mass)
        gas_frac[i] = gas_mass[i] / be_mass
        mean_radius[i] = sum([x_data[i] * mass for i, mass in enumerate(state_dict[state])]) / gas_mass[i]
        #radius_frac[i] = mean_radius[i] / x_data[-1]
        radius_frac[i] = mean_radius[i] / be_radius

    fpath = os.path.join(os.getcwd(), outdir, f"{str(dsfn)}_state_info.h5")
    gas_mass.write_hdf5(fpath, dataset_name='gas_mass')
    gas_frac.write_hdf5(fpath, dataset_name='gas_fraction')
    mean_radius.write_hdf5(fpath, dataset_name='mean_radius')
    radius_frac.write_hdf5(fpath, dataset_name='radius_fraction')




if __name__ == "__main__":

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields= False, load_data=False)
   
    ap.add_operation(create_state_profiles, indir_vel='Profiles/Velocity_and_timescale_profiles', \
                     infile_be= 'Profiles/Bonnor_ebert_profiles/be_time_mass_radius.h5', outdir='Profiles/State_profiles_be')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
