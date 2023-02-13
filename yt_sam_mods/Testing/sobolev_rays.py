import os
import yt
import ytree
yt.enable_parallelism()
import numpy as np

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import *
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods import add_p2p_fields

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc



POSITION = unyt_array([16949.58529203, 18979.38793353, 17409.57691984], 'pc')
OUTFILE = "sobolev.h5"
NUM_RAYS = 20


def ray_trace(node, outfile= "sobolev.h5", num_rays= 20, ray_res= None):

    ds = node.ds
    ds.add_gradient_fields(('gas', 'density'))
    sphere = node.sphere
    halo_center = sphere.center
    star_coords = POSITION
    d_halo = np.linalg.norm((halo_center - star_coords).to('pc')) * pc
    vec_length = d_halo + 2 * sphere.radius
    diff_max = np.arctan(sphere.radius / d_halo)
    x, y, z = halo_center - star_coords
    
    for n in range(num_rays):

        theta = np.arctan(np.sqrt(x**2 + y**2)/z) + np.random.uniform(low= -diff_max, high=diff_max)
        phi = np.arctan(x/y) + np.random.uniform(low= -diff_max, high= diff_max)
        x_ = vec_length * np.cos(phi) * np.cos(theta)
        y_ = vec_length * np.sin(phi) * np.cos(theta)
        z_ = vec_length * np.sin(theta)
        vec_end = transpose_unyt([x_, y_, z_])

        if (ray_res == None):
            ray = ds.r[star_coords:vec_end]
        else :
            #ray = ds.r[star_coords:vec_end:(ray_res)j]
            ray = ds.r[star_coords:vec_end]
        densities = ray[('gas', 'density')][np.argsort(ray['t'])]
        distances = sorted(ray['t']) * vec_length

        grads_xyz = [ray[('gas', f'density_gradient_{ax}')][np.argsort(ray['t'])] for ax in 'xyz']
        grads_vec = transpose_unyt([unyt_quantity(np.dot(ray.vec.to('cm'), transpose_unyt(grad).to('g/cm**4')) \
                                                  / vec_length) for grad in zip(*grads_xyz)])
        grads_vec.write_hdf5("dumbass_wont_let_me_write.h5", dataset_name= "grads")
        l_sob = (max(densities) / grads_vec[np.argmax(densities)]).to('pc')

        l_sob.write_hdf5(outfile, dataset_name= "shock_length")
        (densities / grads_vec).to('pc').write_hdf5(outfile, dataset_name= "sobolev_lengths")
        grads_vec.write_hdf5(outfile, dataset_name= "gradients")
        densities.write_hdf5(outfile, dataset_name= "densities")
        


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
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(ray_trace, outfile= OUTFILE, num_rays= NUM_RAYS)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
