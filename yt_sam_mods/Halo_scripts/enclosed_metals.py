import numpy as np
import os
import yt
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import align_sphere, return_sphere, transpose_unyt
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc




def enclosed_metals(node, outdir='.', nbins=120, take_log= False):

    ds = node.ds
    sphere = node.sphere

    virial_radius = sphere.radius.in_units('pc')
    dx = 2 * ds.index.get_smallest_dx().in_units('pc')
    if (take_log):
        radii = unyt_array(np.logspace(np.log10(dx.value), np.log10(virial_radius.value), num= nbins), 'pc')
    else :
        radii = np.linspace(dx, virial_radius, num= nbins)

    metal_masses = [ds.sphere(sphere.center, radius).quantities.total_quantity(('gas', 'metal3_mass')) for radius in radii]
    metal_masses = transpose_unyt(metal_masses).in_units('Msun')
    metal_ztot = [ds.sphere(sphere.center, radius).quantities.weighted_average_quantity(('gas', 'metallicity3_min7'), ('gas', 'cell_mass')) for radius in radii]
    metal_ztot = transpose_unyt(metal_ztot).in_units('Zsun')

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_enclosed_metals.h5")
    metal_masses.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='metal3_mass')
    metal_ztot.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='metallicity3')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')


def perform_test(node):

    ds = node.ds
    sp = ds.sphere(unyt_array([2000, 2000, 2000], 'pc'), unyt_quantity(200, 'pc'))
    y = sp.quantities.weighted_average_quantity(('gas', 'metallicity3'), ('gas', 'cell_mass'))
    print(y)
    z = sp.quantities.weighted_average_quantity(('gas', 'metallicity3_min7'), ('gas', 'cell_mass'))
    print(z)
    x = sp.quantities.weighted_average_quantity(('gas', 'metal3_mass'), ('gas', 'cell_mass'))
    print(x)



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
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    ap.add_operation(perform_test)

    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)

    ap.add_operation(enclosed_metals, outdir= 'Profiles/Enclosed_metals_lin', nbins= 120)
    ap.add_operation(enclosed_metals, outdir= 'Profiles/Enclosed_metals_log', nbins= 120, take_log= True)

    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
