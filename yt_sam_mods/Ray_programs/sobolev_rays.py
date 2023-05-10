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
from yt.extensions.sam_mods.graph_funcs import get_field_dict
from yt.extensions.sam_mods import add_p2p_fields

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc


#POP3_POSITION = [0.53784, 0.53817, 0.52178] # put these values into graph_funcs.py later!!!!
POP3_POSITION = [0.53805048, 0.53802318, 0.52146589]
OUTDIR = "Sobolev/Ray_profiles"
STRETCH_FACTOR = 1.2
COMPACT_FACTOR = 10.0
NUM_RAYS = 10
RAY_COLOR = 'orange'
BOXSIZE = unyt_quantity(1000.0, 'pc')
FIELDS = [('gas', 'density'), ('gas', 'temperature'), ('gas', 'H2_p0_fraction'), \
          ('gas', 'pressure'), ('gas', 'metallicity3'), ('gas', 'entropy')]
VEC_FIELDS = [('gas', 'velocity')]


def ray_trace(node, outdir= '.', num_rays= 20):

    ray_res= None
    ds = node.ds
    yt.mylog.info(f"Current dataset is {str(ds)}")
    ds.add_gradient_fields(('gas', 'density'))
    sphere = node.sphere
    halo_center = sphere.center
    star_coords = ds.arr(POP3_POSITION, 'code_length')
    d_halo = ds.quan(np.linalg.norm((halo_center - star_coords).to('pc')), 'pc')
    vec_length = d_halo + STRETCH_FACTOR * sphere.radius
    diff_max = np.arctan(sphere.radius / (d_halo * COMPACT_FACTOR))
    x, y, z = halo_center - star_coords

    field_dict = {field[1]: [[] for n in range (num_rays)] for field in FIELDS}
    vec_field_dict = {field[1]: [[] for n in range (num_rays)] for field in VEC_FIELDS}
    distances = [[] for n in range (num_rays)]
    psi = [0 for n in range (num_rays)]
    plots = [yt.ProjectionPlot(ds, ax, ('gas', 'density'), center= halo_center, width= BOXSIZE) for ax in 'xyz']
    
    for n in range(num_rays):
        
        if (n == 0):
            scatter_theta = 0
            scatter_phi = 0
        else :
            scatter_theta = np.random.uniform(low= -diff_max/2, high=diff_max/2)
            scatter_phi = np.random.uniform(low= -diff_max/2, high=diff_max/2)
        theta = np.arctan2(np.sqrt(x**2 + y**2), z) + scatter_theta
        phi = np.arctan2(y, x) + scatter_phi
        x_ = vec_length * np.sin(theta) * np.cos(phi)
        y_ = vec_length * np.sin(theta) * np.sin(phi)
        z_ = vec_length * np.cos(theta)
        vec_end = star_coords + transpose_unyt([x_, y_, z_])
        psi[n] = np.arctan2(np.sqrt(np.tan(phi)**2 + np.tan(theta)**2))
        
        for plt in plots:
            plt.annotate_line(star_coords, vec_end, coord_system='data', plot_args={'color': RAY_COLOR, 'linewidth': 1})
        if (ray_res == None):
            ray = ds.r[star_coords:vec_end]
        else :
            #ray = ds.r[star_coords:vec_end:(ray_res)j]
            ray = ds.r[star_coords:vec_end]

        for field in field_dict.keys():
            field_dict[field][n] = ray[field][np.argsort(ray['t'])]
            yt.mylog.info(f"Calculated ray {n} for {field}")
        for field in vec_field_dict.keys():
            vec_field_dict[field][n] = transpose_unyt([transpose_unyt([ray[f'{field}_{ax}'][np.argsort(ray['t'])][i] \
                                                                       for ax in 'xyz']).dot(ray.vec / vec_length) \
                                                       for i in range (len(ray['t']))])
            yt.mylog.info(f"Calculated ray {n} for {field}")
        distances[n] = sorted(ray['t']) * vec_length

        #grads_xyz = [ray[('gas', f'density_gradient_{ax}')][np.argsort(ray['t'])] for ax in 'xyz']
        #grads_vec = transpose_unyt([unyt_quantity(np.dot(ray.vec / vec_length, transpose_unyt(grad).to('g/cm**4')),
        #                                          'g/cm**4') for grad in zip(*grads_xyz)])
        #sobolev_lengths[n] = (max(densities) / grads_vec[np.argmax(densities)]).to('pc')

    for i, plt in enumerate(plots):
        plt.save(os.path.join(outdir, f"{str(ds)}_projection_{i}_sob_lines.png"))

    outfile = os.path.join(os.getcwd(), outdir, f"{str(ds)}_sobolev.h5")
    field_dict.update(vec_field_dict)
    for field in field_dict.keys():
        for n in range (num_rays):
            field_dict[field][n].write_hdf5(outfile, group_name= f"{field}_{n}")
    for n in range (num_rays):
        distances[n].write_hdf5(outfile, group_name= f"distances_{n}")
    psi.write_hdf5(outfile, group_name="angle_list")



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
    
    ap.add_operation(ray_trace, outdir= OUTDIR, num_rays= NUM_RAYS)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    #N_COMPLETED = 66
    #tree_mod = list(tree['prog'])[N_COMPLETED:]
    #for node in ytree.parallel_trees(tree_mod):
    #    ap.process_target(node)
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
