import statistics
import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.p2p.misc import return_sphere, align_sphere, transpose_unyt
from yt.extensions.p2p.profiles import my_profile
from yt.extensions.p2p.plots import *
from yt.extensions.p2p.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.p2p import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc




def bonnor_ebert(node, outdir='.', nbins=120, take_log= False):

    ds = node.ds
    sphere = node.sphere
    
    dx = ds.index.get_smallest_dx()
    profile = my_profile(sphere, 'radius',
                         [('gas', 'velocity_magnitude'), ('gas', 'density'), ('gas', 'cooling_time'), ('gas', 'pressure')],
                         logs={'radius': take_log}, n_bins= nbins, weight_field= ('gas', 'cell_mass'))
    radii = profile.x_bins[1:].in_units('pc')
    spheres = [ds.sphere(sphere.center, radius) for radius in radii]
    mtot_prof = [sum(sphere.quantities.total_mass()) for sphere in spheres]
    #mtot_prof = [sphere.quantities.total_mass()[0] for sphere in spheres]
    cs_prof = [sphere.quantities.weighted_average_quantity(('gas', 'sound_speed'), ('gas', 'cell_mass')) for sphere in spheres]
    vsigma_prof = [sphere.quantities.weighted_average_quantity(('gas', 'vsigma'), ('gas', 'cell_mass')) for sphere in spheres]
    
    mtot_prof = transpose_unyt(mtot_prof).in_units('Msun')
    cs_prof = transpose_unyt(cs_prof).in_units('km/s')
    vsigma_prof = transpose_unyt(vsigma_prof).in_units('km/s')
    
    tff_norm_ratio = ( (np.pi / np.sqrt(8 * G)) * cs_prof * np.sqrt(radii / mtot_prof) ).in_units('')
    tff_turb_ratio = ( (np.pi / np.sqrt(8 * G)) * vsigma_prof * np.sqrt(radii / mtot_prof) ).in_units('')
    tcool_norm_ratio = ( cs_prof * profile[('gas', 'cooling_time')] / radii ).in_units('')
    tcool_turb_ratio = ( vsigma_prof * profile[('gas', 'cooling_time')] / radii ).in_units('')
    #be_prof = mtot_prof - ( (1.18 * cs_prof**4) / (profile[('gas', 'velocity_magnitude')] * np.sqrt(2 * profile[('gas', 'density')] * G**3)) )
    be_prof = mtot_prof - ( (1.18 * cs_prof**4) / np.sqrt(G**3 * profile[('gas', 'pressure')]) ) 
    
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_bonnor_ebert.h5")
    be_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='deviation')
    tff_norm_ratio.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='tff_ratio_norm')
    tff_turb_ratio.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='tff_ratio_turb')
    tcool_norm_ratio.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='tcool_ratio_norm')
    tcool_turb_ratio.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='tcool_ratio_turb')
    mtot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='mass')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')
    



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
    ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(bonnor_ebert, outdir= 'Profiles/Bonnor_ebert_lin', nbins= 120)
    ap.add_operation(bonnor_ebert, outdir= 'Profiles/Bonnor_ebert_log', nbins= 120, take_log=True)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
