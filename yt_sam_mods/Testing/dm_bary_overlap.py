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
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc



def density_overlap(node, outdir='.', nbins=120):

    ds = node.ds
    sphere = node.sphere

    #ds.add_gradient_fields(('gas', 'density'))
    #ds.add_gradient_fields(('gas', 'dark_matter_density'))
    #profile_test = my_profile(sphere, ('index', 'radius'), [('gas', 'vortical_time')], logs={'radius': True}, n_bins= nbins, weight_field= ('gas', 'cell_mass'))
    profile = my_profile(sphere, ('index', 'radius'), [('gas', 'overlap')], logs={'radius': True}, n_bins= nbins, weight_field= ('gas', 'cell_mass'))
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_entropy_profile.h5")
    profile.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='overlap')





if __name__ == "__main__":

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    #tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    #ap.add_operation(yt_dataset, data_dir, es)
    #ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(density_overlap, outdir= 'Profiles/Overlap_profiles')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    
    tree = a[0]

    '''
    N_MISSED = 19
    tree_mod = list(tree['prog'])[N_MISSED:]
    for node in ytree.parallel_trees(tree_mod):
        ap.process_target(node)
    '''

    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
