import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import return_sphere_pop3, align_sphere, modify_grackle
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods.graph_funcs import get_field_dict
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc



BIN_DENSITY = 20
BOXSIZE = unyt_quantity(300, 'pc')

FIELDS = ['density', 'H2_p0_fraction', 'metallicity3', 'temperature']


if __name__ == "__main__":

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    #tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields= True)
    #ap.add_operation(yt_dataset, data_dir, es)
    #ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere_pop3, 'pisn')

    field_info = dict([(('gas', field), (get_field_dict(field)['colormap'], \
                                         get_field_dict(field)['units'], \
                                         get_field_dict(field)['limits'])) for field in FIELDS])
    #ap.add_operation(region_projections, 'pisn', fields, size= BOXSIZE, fixed_size=True, output_dir='Projections_pop3')
    ap.add_operation(region_projections, 'pisn', field_info, output_dir='Projections_pop3')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    #N_MISSED = 23
    #tree_mod = list(tree['prog'])[N_MISSED:]
    #for node in ytree.parallel_trees(tree_mod):
    #    ap.process_target(node)
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
