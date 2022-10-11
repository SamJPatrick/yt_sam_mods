import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import return_sphere, align_sphere, transpose_unyt
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc




def bonnor_ebert(node, veldir= '.', normdir= '.', outdir='.'):

    dsfn = node.ds_filename.split('/')[-1]
    normdf_w = yt.load(os.path.join(normdir, f"{str(dsfn)}_profile_weight_field_cell_mass.h5"))
    normdf_u = yt.load(os.path.join(normdir, f"{str(dsfn)}_profile_weight_field_None.h5"))
    veldf = yt.load(os.path.join(veldir, f"{str(dsfn)}_1D_profile_radius_cell_mass.h5"))
    
    used = normdf_u.profile.used
    x_data = normdf_u.profile.x[used].in_units('pc')
    mass = normdf_u.profile[('data', 'matter_mass')][used].in_units('Msun')
    tot_mass = [sum(mass[:i+1]) for i in range (len(x_data))]
    tot_mass = transpose_unyt(tot_mass)
    gas_mass = normdf_u.profile[('data', 'cell_mass')][used].in_units('Msun')
    tot_gas_mass = [sum(gas_mass[:i+1]) for i in range (len(x_data))]
    tot_gas_mass = transpose_unyt(tot_gas_mass)

    prefac = (np.pi / np.sqrt(8 * G))
    cs = veldf.profile[('data', 'sound_speed')][used]
    vel = veldf.profile[('data', 'velocity_magnitude')][used]
    sigma = veldf.profile.standard_deviation[('data', 'velocity_magnitude')][used]
    turb_cs = sigma
    #turb_cs = np.sqrt(sigma**2 + cs**2)
    tff_norm_ratio = (prefac * cs * np.sqrt(x_data / tot_mass)).in_units('')
    tff_turb_ratio = (prefac * turb_cs * np.sqrt(x_data / tot_mass)).in_units('')
    tcool_norm_ratio = (cs * veldf.profile[('data', 'cooling_time')][used] / x_data).in_units('')
    tcool_turb_ratio = (turb_cs * veldf.profile[('data', 'cooling_time')][used] / x_data).in_units('')
    pressure_norm = pressure_turb = normdf_w.profile[('data', 'pressure')][used]
    #pressure_norm = normdf_w.profile[('data', 'pressure')][used] + vel**2 * normdf_w.profile[('data', 'density')][used]
    #pressure_turb = normdf_w.profile[('data', 'pressure')][used] + turb_cs**2 * normdf_w.profile[('data', 'density')][used]
    be_norm = (1.18 * cs**4) / np.sqrt(G**3 * pressure_norm)
    be_turb = (1.18 * turb_cs**4) / np.sqrt(G**3 * pressure_turb)
    be_norm_prof = tot_mass / be_norm
    be_turb_prof = tot_mass / be_turb
    
    fpath = os.path.join(os.getcwd(), outdir, f"{str(dsfn)}_bonnor_ebert.h5")
    be_norm_prof.write_hdf5(fpath, group_name='be_ratio_norm')
    be_turb_prof.write_hdf5(fpath, group_name='be_ratio_turb')
    tff_norm_ratio.write_hdf5(fpath, group_name='tff_ratio_norm')
    tff_turb_ratio.write_hdf5(fpath, group_name='tff_ratio_turb')
    tcool_norm_ratio.write_hdf5(fpath, group_name='tcool_ratio_norm')
    tcool_turb_ratio.write_hdf5(fpath, group_name='tcool_ratio_turb')
    tot_mass.write_hdf5(fpath, group_name='total_mass')
    tot_gas_mass.write_hdf5(fpath, group_name='total_gas_mass')
    x_data.write_hdf5(fpath, group_name='radius')





if __name__ == "__main__":

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")
    
    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=False, load_data=False)    
    ap.add_operation(bonnor_ebert, veldir= 'Profiles/Velocity_and_timescale_profiles', \
                     normdir= 'Profiles/Normal_profiles', outdir= 'Profiles/Bonnor_ebert_profiles')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
