import sys
import os
import yt
import h5py
yt.enable_parallelism()
import ytree
import numpy as np

from ytree.analysis import AnalysisPipeline
from yt.extensions.sam_mods.misc import return_sphere, align_sphere
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods import add_p2p_fields


from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc


PRE_DIR = "Britton_sim_data"
VELDIR = "Profiles/Velocity_profiles"
NORMDIR = "Profiles/Normal_profiles"
OUTDIR = "Profiles/Pressure_profiles"


def pressures(node, veldir= '.', normdir= '.', outdir='.'):

    dsfn = node.ds_filename.split('/')[-1]
    normdf_w = yt.load(os.path.join(normdir, f"{str(dsfn)}_profile_weight_field_cell_mass.h5"))
    normdf_u = yt.load(os.path.join(normdir, f"{str(dsfn)}_profile_weight_field_None.h5"))
    veldf = yt.load(os.path.join(veldir, f"{str(dsfn)}_profile_weight_field_cell_mass.h5"))
    
    used = veldf.profile.used
    x_data = veldf.profile.x[used].to("pc")
    mass = normdf_u.profile[('data', 'matter_mass')][used].in_units('Msun')
    tot_mass = [sum(mass[:i+1]) for i in range (len(x_data))]
    
    p_hydro_dx = G * normdf_w.profile[('data', 'density')][used] * tot_mass / x_data
    p_hydro = [sum(p_hydro_dx[:i+1]) for i in range(len(x_data), 0, -1)]
    p_hydro = transpose_unyt(p_hydro)
    p_ram = normdf_w.profile[('data', 'density')][used] * veldf.profile[('data', 'velocity_spherical_radius')][used]**2
    p_turb = normdf_w.profile[('data', 'density')][used] * veldf.profile.standard_deviation[('data', 'velocity_magnitude')][used]**2
    p_therm = normdf_w.profile[('data', 'pressure')][used]
    
    fpath = os.path.join(outdir, f"{str(dsfn)}_pressures.h5")
    p_hydro.in_units('Pa').write_hdf5(fpath, group_name='hydrostatic')
    p_ram.in_units('Pa').write_hdf5(fpath, group_name='ram')
    p_turb.in_units('Pa').write_hdf5(fpath, group_name='turbulent')
    p_therm.in_units('Pa').write_hdf5(fpath, group_name='thermal')
    mass.write_hdf5(fpath, group_name='mass')
    x_data.write_hdf5(fpath, group_name='radius')





if __name__ == "__main__":

    try:
        star_mode = sys.argv[1].upper()
    except IndexError:
        print("Error, 'star_mode' argument required")
        quit()

    if (star_mode == 'PISN'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"
    elif (star_mode == 'HN'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"
    elif (star_mode == 'CCSN'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"
    else :
        print("Error, 'star_mode' not recognised")
        quit()
    veldir = os.path.join(PRE_DIR, star_mode, VELDIR)
    normdir = os.path.join(PRE_DIR, star_mode, NORMDIR)
    outdir = os.path.join(PRE_DIR, star_mode, OUTDIR)

    #es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")
    
    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=False, load_data=False)   
    ap.add_operation(pressures, veldir= veldir, normdir= normdir, outdir= outdir) 
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
