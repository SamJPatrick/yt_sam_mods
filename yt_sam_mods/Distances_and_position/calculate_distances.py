import numpy as np
import os
import yt
yt.enable_parallelism()
import ytree

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import retrieve_positions
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh




def print_distances(node):

    star_position = node.star_position
    halo_position = node.halo_position
    distance = node.distance
    yt.mylog.info(f"Star position: {star_position}")
    yt.mylog.info(f"Halo position: {halo_position}")
    yt.mylog.info(f"Distance: {distance}")




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
    ap.add_operation(yt_dataset, data_dir, add_fields= True)
    #ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(retrieve_positions, position_file= "positions_info.h5")
    ap.add_operation(print_distances)

    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
