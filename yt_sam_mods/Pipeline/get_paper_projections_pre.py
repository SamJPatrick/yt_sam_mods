import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

#from yt.extensions.sam_mods.misc import sphere_icom
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods.graph_funcs import get_field_dict
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods import add_p2p_fields
from yt.utilities.math_utils import ortho_find

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc


SIM_PATH = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large"
TREE_PATH = "merger_trees/target_halos/target_halos.h5"
OUTDIR =  '.'

POP3_POSITION = [0.53805048, 0.53802318, 0.52146589]
BOXSIZE_PROJ = unyt_quantity(300, 'pc')
BOXSIZE_SLICE = unyt_quantity(300, 'pc')

FIELDS_SLICE = ['density', 'temperature']
FIELDS_PROJ = ['density', 'temperature']
limit_dict = {'density': unyt_array([1e-27, 1e-21], 'g/cm**3'),
              'temperature': unyt_array([1e1, 1e8], 'K'),
              'metallicity3': unyt_array([1e-7, 1e0], 'Zsun')}
cmap_dict = {'density': 'turbo', 'temperature': 'inferno', 'metallicity3': 'cool'}
big_dict = {'ccsn': {'sim_dir': "cc_512_collapse_solar_dust",
                     'dump_nums_slice': [104],
                     'dump_nums_proj': []},
            'hn' : {'sim_dir': "hyper_512_collapse_solar_dust",
                    'dump_nums_slice': [],
                    'dump_nums_proj': []},
            'pisn': {'sim_dir': "pisn_solo",
                     'dump_nums_slice': [130],
                     'dump_nums_proj': []}}




def get_node(tree, time):

    time_diffs = transpose_unyt([np.abs(node['time'] - time) for node in tree['prog']])
    node = list(tree['prog'])[np.argmin(time_diffs)]
    return node


def sphere_icom(ds, sphere, inner_radius, stepsize=0.5, com_kwargs=None):
    
    if com_kwargs is None:
        com_kwargs = {}

    yield sphere
    while (sphere.radius > inner_radius):
        com = sphere.quantities.center_of_mass(**com_kwargs)
        yt.mylog.info(f"Current center of mass is {com.to('pc')}")
        try:
            sphere = ds.sphere(com, (1 - stepsize) * sphere.radius)
            yield sphere
        except YTSphereTooSmall:
            yield None
            break


def get_sphere(node, ds):
    
    radius = ds.quan(node["virial_radius"], "unitary").to('pc')
    #if "icom_gas2_position_x" in node.arbor.field_list:
    #    center = ds.arr(node["icom_gas2_position"], "unitary").to('pc')
    #    sphere = ds.sphere(radius, center)
    #else:
    center = ds.arr(node["position"], "unitary").to('pc')
    sphere = ds.sphere(center, radius)
    for new_sphere in sphere_icom(ds, sphere, 4*sphere[("gas", "dx")].min(),
                                  com_kwargs={'use_particles':False, 'use_gas':True}):
        if new_sphere is None:
            break
    return new_sphere


def get_ray(ds, halo_center):

    star_coords = ds.arr(POP3_POSITION, 'code_length').to('pc')
    d_halo = ds.quan(np.linalg.norm((halo_center - star_coords).to('pc')), 'pc')
    unit_vec = (star_coords - halo_center) / d_halo
    ray = ds.r[star_coords:halo_center]
    yt.mylog.info(f"Inter-halo distance for dataset {str(ds)} is {d_halo}")
    ray_vec = ds.arr(ray.vec, 'code_length').to('pc')
    north_vecs = [ds.arr(vec, 'code_length').to('pc') for vec in ortho_find(ray_vec)[1:]]
    ray_center = ray.start_point.to('pc') + ray_vec / 2
    return ray_center, ray_vec, north_vecs





if __name__ == "__main__":
    
    for star_type in list(big_dict.keys()):
        tree_path = os.path.join(SIM_PATH, big_dict[star_type]['sim_dir'], TREE_PATH)
        a = ytree.load(tree_path)

        dir_paths = [os.path.join(SIM_PATH, big_dict[star_type]['sim_dir'], f"DD{dump_num:04d}/DD{dump_num:04d}")
                     for dump_num in big_dict[star_type]['dump_nums_proj']]
        for dir_path in dir_paths:
            ds = yt.load(dir_path)
            add_p2p_fields(ds)
            node = get_node(a[0], ds.current_time)
            #sphere = get_sphere(node, ds)
            sphere = ds.sphere(ds.arr(node['position'], 'unitary').to('pc'), ds.quan(node['virial_radius'], 'unitary').to('pc'))
            halo_center = sphere.center
            for field in FIELDS_PROJ:
                p = yt.ProjectionPlot(ds, "x", [("gas", field)], center= halo_center, width= BOXSIZE_PROJ, weight_field=("gas", "density"))
                p.set_cmap(field, cmap_dict[field])
                p.set_axes_unit('pc')
                p.set_zlim(field, *limit_dict[field])
                p.frb.save_as_dataset(filename= os.path.join(OUTDIR, f"{str(ds)}_{field}_proj.h5"))

        dir_paths = [os.path.join(SIM_PATH, big_dict[star_type]['sim_dir'], f"DD{dump_num:04d}/DD{dump_num:04d}")
                     for dump_num in big_dict[star_type]['dump_nums_slice']]
        for dir_path in dir_paths:
            ds = yt.load(dir_path)
            add_p2p_fields(ds)
            node = get_node(a[0], ds.current_time)
            #sphere = get_sphere(node, ds)
            sphere = ds.sphere(ds.arr(node['position'], 'unitary').to('pc'), ds.quan(node['virial_radius'], 'unitary').to('pc'))
            halo_center = sphere.center
            #ray_center, ray_vec, north_vecs = get_ray(ds, halo_center)
            for field in FIELDS_SLICE:
                p = yt.SlicePlot(ds, "x", [("gas", field)], center= halo_center, width= BOXSIZE_PROJ)
                #p = yt.OffAxisSlicePlot(ds, north_vecs[1], [('gas', field)], center= ray_center,
                #                        width= ds.quan(BOXSIZE_SLICE.value, BOXSIZE_SLICE.units), north_vector= ray_vec)
                p.set_cmap(field, cmap_dict[field])
                p.set_axes_unit('pc')
                p.set_zlim(field, *limit_dict[field])
                p.frb.save_as_dataset(filename= os.path.join(OUTDIR, f"{str(ds)}_{field}_slice.h5"))
