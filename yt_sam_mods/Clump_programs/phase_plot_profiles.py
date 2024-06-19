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
OUTDIR = "Profiles"
UNIT_DICT = {('gas', 'number_density'): 'cm**(-3)', ('gas', 'temperature'): 'K', ('index', 'radius'): 'pc'}


def mass_weighted_profiles(node, outdir="."):

    #bin_fields = [('gas', 'number_density'), ('gas', 'temperature')]
    bin_fields = [('gas', 'number_density'), ('gas', 'temperature'), ('index', 'radius')]
    profile_fields = [('gas', 'cell_mass')]
    weight_field = None
    profile = my_profile(node.sphere, bin_fields, profile_fields, units= UNIT_DICT,
                         weight_field=weight_field, accumulation=False, bin_density= BIN_DENSITY)
    fn = os.path.join(outdir, f"{str(node.ds)}_profile_number_density_temperature_radius.h5")
    profile.save_as_dataset(filename=fn)





if __name__ == "__main__":

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    #tree_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

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
    ap.process_target(tree)
