import os
import sys
import yt
import ytree
yt.enable_parallelism()
import numpy as np

from ytree.analysis import AnalysisPipeline
from yt.utilities.math_utils import ortho_find

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
OUTDIR_PROF = "Sobolev/Ray_profiles"
OUTDIR_PLOT = "Sobolev/Ray_profiles"
OUTDIR_PROJ = "Slices_ray"

STRETCH_FACTOR = 1.2
COMPACT_FACTOR = 10.0
NUM_RAYS = 10
RAY_COLOR = 'orange'
BOXSIZE = (1000.0, 'pc')

FIELDS = ['density', 'dark_matter_density', 'H2_p0_fraction', 'El_fraction', 'cooling_time', \
          'temperature', 'pressure', 'metallicity3', 'entropy', 'cell_mass']
PLOT_FIELDS = ['density', 'El_fraction', 'cooling_time', \
               'temperature', 'pressure', 'cooling_intensity', 'cell_mass']
VEC_FIELDS = ['velocity']




def pack_fields(distance_arr, field_dict):

    num_rays = len(distance_arr)
    for field, field_arr in field_dict.items():
        assert len(field_arr) == num_rays, f"Error, {field} does not contain the appropriate {num_rays} of rays"
    min_length = np.min([distances[-1].to('pc') for distances in distance_arr])
    max_bins = np.max([len(distances) for distances in distance_arr])
    ref_ray = unyt_array(np.linspace(0.0, min_length, num= max_bins), 'pc')
    
    packed_dict = {field: [unyt_array(np.zeros(max_bins), get_field_dict(field)['units']) for num in range (num_rays)] \
                   for field in field_dict.keys()}
    for field in field_dict.keys():
        for num in range (num_rays):
            for i in range (max_bins):
                if (ref_ray[i] <= distance_arr[num][0]):
                    ind1 = 0
                else :
                    ind1 = np.argwhere((ref_ray[i] - distance_arr[num]) > 0)[-1].item()
                if (ref_ray[i] >= distance_arr[num][-1]):
                    ind2 = len(distance_arr[num]) - 1
                else :
                    ind2 = np.argwhere(distance_arr[num] - ref_ray[i] > 0)[0].item()
                packed_dict[field][num][i] = field_dict[field][num][ind1] + \
                    ((field_dict[field][num][ind2] - field_dict[field][num][ind1]) \
                     / (distance_arr[num][ind2] - distance_arr[num][ind1])) * (ref_ray[i] -  distance_arr[num][ind1])

    return ref_ray, packed_dict




def get_rays(node, num_rays= 10):

    ds = node.ds
    sphere = node.sphere
    halo_center = sphere.center.to('pc')
    star_coords = ds.arr(POP3_POSITION, 'code_length').to('pc')
    d_halo = ds.quan(np.linalg.norm((halo_center - star_coords).to('pc')), 'pc')
    
    unit_vec = (star_coords - halo_center) / d_halo
    master_ray = ds.r[star_coords:halo_center]
    distances = sorted(master_ray['t']) * d_halo
    masses = master_ray['cell_mass'][np.argsort(master_ray['t'])]
    ave_distance = np.average(distances, weights= masses)
    ray_com = star_coords + unit_vec * ave_distance
    yt.mylog.info(f"Inter-halo distance for dataset {str(ds)} is {d_halo}")
    yt.mylog.info(f"COM distance along ray for dataset {str(ds)} is {ave_distance} from PopIII star at {ray_com}")

    diff_max = np.arctan(sphere.radius / (d_halo * COMPACT_FACTOR))
    vec_length = d_halo + STRETCH_FACTOR * sphere.radius
    x, y, z = halo_center - star_coords
    probe_psis = unyt_array(np.zeros(num_rays), 'dimensionless')
    probe_rays = [None for n in range (num_rays)]
    for n in range(num_rays):
        if (n == 0):
            scatter_theta = 0
            scatter_phi = 0
        else :
            scatter_theta = np.random.uniform(low= -diff_max/2, high=diff_max/2)
            scatter_phi = np.random.uniform(low= -diff_max/2, high=diff_max/2)
        theta = np.arctan2(np.sqrt(x**2 + y**2), z) + scatter_theta
        phi = np.arctan2(y, x) + scatter_phi
        yt.mylog.info(f"Angles: {theta} & {phi}")
        x_ = vec_length * np.sin(theta) * np.cos(phi)
        y_ = vec_length * np.sin(theta) * np.sin(phi)
        z_ = vec_length * np.cos(theta)
        vec_end = star_coords + transpose_unyt([x_, y_, z_])
        probe_rays[n] = ds.r[star_coords:vec_end]
        probe_psis[n] = np.arctan(np.sqrt(np.tan(phi)**2 + np.tan(theta)**2))

    node.master_ray = master_ray
    node.probe_psis = probe_psis
    node.probe_rays = probe_rays




def profile_rays(node, outdir= '.'):

    rays = node.probe_rays
    num_rays = len(rays)
    field_dict = {field: [[] for n in range (num_rays)] for field in FIELDS}
    vec_field_dict = {**{f'{field}_para': [[] for n in range (num_rays)] for field in VEC_FIELDS}, \
                      **{f'{field}_norm': [[] for n in range (num_rays)] for field in VEC_FIELDS}}
    distances = [[] for n in range (num_rays)]
    
    for n, ray in enumerate(rays):
        ortho_vec = ortho_find(ray.vec.to('pc'))[1]
        yt.mylog.info(f"Ray vector is {(ray.vec).to('pc')}")
        yt.mylog.info(f"Orthogonal vector is {ortho_vec}")
        ray_length = unyt_quantity(np.linalg.norm((ray.vec).to('pc')), 'pc')
        for field in field_dict.keys():
            field_dict[field][n] = ray[field][np.argsort(ray['t'])].to(get_field_dict(field)['units'])
            yt.mylog.info(f"Calculated ray {n} for {field}")
        for field in VEC_FIELDS:
            vec_field_dict[f'{field}_para'][n] = transpose_unyt([transpose_unyt([ray[f'{field}_{ax}'][np.argsort(ray['t'])][i] \
                                                                                 for ax in 'xyz']).dot(ray.vec / ray_length) \
                                                                 for i in range (len(ray['t']))])
            vec_field_dict[f'{field}_norm'][n] = transpose_unyt([transpose_unyt([ray[f'{field}_{ax}'][np.argsort(ray['t'])][i] \
                                                                                 for ax in 'xyz']).dot(ortho_vec) \
                                                                 for i in range (len(ray['t']))])
            yt.mylog.info(f"Calculated ray {n} for {field}")
        distances[n] = sorted(ray['t']) * ray_length
    field_dict.update(vec_field_dict)

    outfile = os.path.join(os.getcwd(), outdir, f"{str(node.ds)}_sobolev.h5")
    for field in field_dict.keys():
        for n in range (num_rays):
            field_dict[field][n].write_hdf5(outfile, group_name= f"{field}_{n}")
    for n in range (num_rays):
        distances[n].write_hdf5(outfile, group_name= f"distances_{n}")
    node.probe_psis.write_hdf5(outfile, group_name="angle_list")

    distances_packed, packed_dict = pack_fields(distances, field_dict)
    outfile = os.path.join(os.getcwd(), outdir, f"{str(node.ds)}_packed.h5")
    for field in packed_dict.keys():
        for n in range (num_rays):
            packed_dict[field][n].write_hdf5(outfile, group_name= f"{field}_{n}")
    distances_packed.write_hdf5(outfile, group_name= f"distances")
    node.probe_psis.write_hdf5(outfile, group_name="angle_list")




def plot_rays(node, star_mode, outdir= '.', color= 'orange', width= 1):

    rays = node.probe_rays
    ds = node.ds
    halo_center = node.sphere.center.to('pc')
    plots = [yt.ProjectionPlot(ds, ax, ('gas', 'density'), center= halo_center, width= BOXSIZE) for ax in 'xyz']
    field_dict = get_field_dict(('gas', 'density'))
    for i, plt in enumerate(plots):
        for ray in rays:
            plt.annotate_line(ray.start_point, ray.end_point, coord_system='data', plot_args={'color': color, 'linewidth': width})
        plt.set_cmap(('gas', 'density'), 'viridis')
        plt.set_unit(('gas', 'density'), 'g/cm**2')
        plt.set_axes_unit('pc')
        plt.set_zlim(('gas', 'density'), 1e-4, 1e-1)
        decorate_plot(node, plt, star_mode)
        plt.save(os.path.join(outdir, f"{str(ds)}_projection_{i}_sob_lines.png"))




def ray_projections(node, star_mode, outdir= '.', num_slices= 20):

    ray = node.master_ray
    ds = node.ds
    ray_vec = ds.arr(ray.vec, 'code_length').to('pc')
    north_vecs = ortho_find(ray.vec)[1:]
    
    if ((type(num_slices) is int) and num_slices > 0):
        centers_para = [ray.start_point + ray.vec * dx for dx in np.linspace(0.0, 1.0, num= num_slices)]
    else :
        centers_para = []
    center_norm = ray.start_point + ray.vec / 2
    yt.mylog.info(f"Centers for dataset {str(ds)} at {ds.current_time} are {centers_para} with center at {center_norm}")

    outdir_para = outdir + "_para"
    outdir_norm = outdir + "_norm"
    for field in PLOT_FIELDS:
        field_dict = get_field_dict(field)
        plots = []
        for center in centers_para:
            plots.append(yt.OffAxisSlicePlot(ds, ray.vec, ('gas', field), center= center, width= BOXSIZE, north_vector= north_vecs[1]))
        plots.append(yt.OffAxisSlicePlot(ds, north_vecs[1], ('gas', field), center= center_norm, width= BOXSIZE, north_vector= ray_vec))
        for i, plt in enumerate(plots):
            plt.set_cmap(field, field_dict['colormap'])
            plt.set_unit(field, field_dict['units'])
            plt.set_axes_unit('pc')
            plt.set_zlim(field, *field_dict['limits'])
            decorate_plot(node, plt, star_mode)
            if (i == len(plots) - 1):
                plt.save(os.path.join(outdir_para, f"{str(ds)}_offax-slice_para_{field}.png"))
            else :
                plt.save(os.path.join(outdir_norm, f"{str(ds)}_offax-slice_{i}_{field}.png"))





if __name__ == "__main__":

    try :
        star_type = sys.argv[1]
    except IndexError:
        star_type = "pisn"
        pass

    if (star_type == 'ccsn'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"
    elif (star_type == 'pisn'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"
    elif (star_type == 'hn'):
        output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/minihalo_analysis"
        sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/simulation.h5"
        data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust"
        tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    #es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    ap.add_operation(get_rays, num_rays= NUM_RAYS)
    
    ap.add_operation(profile_rays, outdir= OUTDIR_PROF)
    ap.add_operation(plot_rays, star_type, outdir= OUTDIR_PLOT, color= RAY_COLOR)
    ap.add_operation(ray_projections, star_type, outdir= OUTDIR_PROJ, num_slices= None)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    #N_COMPLETED = 11
    #tree_mod = list(tree['prog'])[N_COMPLETED:]
    #for node in ytree.parallel_trees(tree_mod):
    #    ap.process_target(node)
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
