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
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, region_projections, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc

from grid_figure import GridFigure



BIN_DENSITY = 20

DENSITY_LIMS = (1e-27, 1e-21)
TEMPERATURE_LIMS = (1e1, 1e5)
METALLICITY_LIMS = (1e-9, 1e-1)
H2_FRACTION_LIMS = (1e-4, 5e-3)




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
    ap.add_operation(yt_dataset, data_dir, add_fields= True)
    #ap.add_operation(yt_dataset, data_dir, es)
    #ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)

    fields = {('gas', 'density'): ('viridis', 'g/cm**3', DENSITY_LIMS),
              ('gas', 'H2_p0_fraction'): ('cool', '', H2_FRACTION_LIMS),
              ('gas', 'metallicity3'): ('cool', '', METALLICITY_LIMS),
              ('gas', 'temperature'): ('plasma', 'K', TEMPERATURE_LIMS)}
    ap.add_operation(region_projections, fields, fixed_size=True, output_dir='Projections')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
