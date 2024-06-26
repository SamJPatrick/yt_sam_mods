import os
import yt
import ytree
yt.enable_parallelism()
import numpy as np

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import *
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, \
    node_profile, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc



BIN_DENSITY = 20
OUTDIR = "Profiles/Normal_profiles"



def sphere_radial_profiles(node, fields, weight_field=None, outdir="."):

    if weight_field is None:
        weight_name = "None"
    else:
        if not (isinstance(weight_field, tuple) and len(weight_field) == 2):
            raise ValueError("weight_field must be a tuple of length 2.")
        weight_name = weight_field[1]

    fn = os.path.join(outdir, f"{str(node.ds)}_profile_weight_field_{weight_name}.h5")
    if os.path.exists(fn):
        return

    data_source = node.sphere
    profile = my_profile(data_source, ("index", "radius"), fields, units={("index", "radius"): "pc"}, \
                         weight_field=weight_field, accumulation=False, bin_density= BIN_DENSITY)
    profile.save_as_dataset(filename=fn)


def mass_weighted_profiles(node, outdir="."):

    fields_weighted = [('gas', 'density'),
                       ('gas', 'dark_matter_density'),
                       ('gas', 'matter_density'),
                       ('gas', 'metallicity3'),
                       ('gas', 'H2_p0_fraction'),
                       ('gas', 'temperature'),
                       ('gas', 'pressure'),
                       ('gas', 'entropy'),
                       ('gas', 'sound_speed'),
                       ('gas', 'accretion_rate'),
                       ('gas', 'accretion_rate_z'),
                       ('gas', 'sound_speed'),
                       ('gas', 'vortical_time'),
                       ('gas', 'cooling_time'),
                       ('gas', 'dynamical_time'),
                       ('gas', 'velocity_magnitude'),
                       ('gas', 'total_dynamical_time')]
    fields_raw = [('gas', 'metal3_mass'),
                  ('gas', 'cell_mass'),
                  ('gas', 'matter_mass'),
                  ('gas', 'dark_matter_mass')]
    sphere_radial_profiles(node, fields_weighted, weight_field=("gas", "cell_mass"), outdir= outdir)
    sphere_radial_profiles(node, fields_raw, weight_field=None, outdir= outdir)






if __name__ == "__main__":

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    tree_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    #tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    #es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    #ap.add_operation(yt_dataset, data_dir, es)
    #ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(mass_weighted_profiles, outdir= OUTDIR)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    # N_MISSED = 19
    # tree_mod = list(tree['prog'])[N_MISSED:]
    # for node in ytree.parallel_trees(tree_mod):
    #     ap.process_target(node)
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
