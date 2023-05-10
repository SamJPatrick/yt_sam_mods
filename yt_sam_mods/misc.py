import yt
import pygrackle
import os
import h5py
import numpy as np
from unyt import unyt_array, unyt_quantity
from unyt import pc, Myr
from yt.extensions.astro_analysis.halo_analysis.halo_catalog.halo_callbacks import periodic_distance
from yt.utilities.exceptions import YTSphereTooSmall
from yt.extensions.sam_mods.graph_funcs import get_sn_energy, get_time_offset



SAFETY_FACTOR = 1.2
POP3_POSITION = [0.53784, 0.53817, 0.52178] # put these values into graph_funcs.py later!!!!
RHO_BKGRD = unyt_quantity(1e-25, 'g/cm**3')
BOXSIZE_DEFAULT = unyt_quantity(300.0, 'pc')
BOXSIZE_MIN = unyt_quantity(100.0, 'pc')



def dirname(path, up=0):
    return "/".join(path.split('/')[:-up-1])


def iterate_center_of_mass(sphere, inner_radius, stepsize=0.05,
                           com_kwargs=None):
    """
    Return a series of spheres centered on the current center of
    mass with radius decreasing by stepsize. Stop when inner_radius
    is reached.
    """

    if com_kwargs is None:
        com_kwargs = {}

    yield sphere
    while (sphere.radius > inner_radius):
        com = sphere.quantities.center_of_mass(**com_kwargs)
        try:
            sphere = sphere.ds.sphere(com, (1-stepsize) * sphere.radius)
            yield sphere
        except YTSphereTooSmall:
            yield None
            break


def sphere_icom(sphere, inner_radius, stepsize=0.1,
                com_kwargs=None, verbose=True):
    center_orig = sphere.center
    old_center = center_orig

    for new_sphere in iterate_center_of_mass(
            sphere, inner_radius, stepsize, com_kwargs):
        if new_sphere is None:
            break

        new_center = new_sphere.center
        if verbose:
            diff = uperiodic_distance(
                old_center.to("unitary"),
                new_center.to("unitary"))
            yt.mylog.info(
                    "Radius: %s, center: %s, diff: %s." %
                    (new_sphere.radius.to('pc'),
                     new_sphere.center.to('unitary'),
                     diff.to('pc')))

        sphere = new_sphere
        old_center = new_center

    if verbose:
        distance = uperiodic_distance(
            center_orig.to("unitary"),
            new_center.to("unitary"))
        yt.mylog.info("Recentering sphere %s away." %
                       distance.to("pc"))

    return sphere.center


def return_sphere(node):
    ds = node.ds
    radius = reunit(ds, node["virial_radius"], "unitary")
    if "icom_gas2_position_x" in node.arbor.field_list:
        center = reunit(ds, node["icom_gas2_position"], "unitary")
    else:
        center = reunit(ds, node["position"], "unitary")
        sphere = ds.sphere(center, radius)
        center = sphere_icom(sphere, 4*sphere["gas", "dx"].min(), \
                             com_kwargs=dict(use_particles=False, use_gas=True))
    #yt.mylog.info(f"TARGET HALO CENTER: {center}")
    sphere = ds.sphere(center, radius)
    node.sphere = sphere


def return_sphere_pop3(node, star_mode):
    ds = node.ds
    center = ds.arr(POP3_POSITION, 'code_length')
    radius_min = BOXSIZE_MIN * (1.0/1.05)
    if (ds.current_time < get_time_offset(star_mode)):
        #radius = reunit(ds, node["virial_radius"], "unitary")
        radius = radius_min
    else :
        st_radius = SAFETY_FACTOR * (get_sn_energy(star_mode) / RHO_BKGRD)**(1/5) *\
            (ds.current_time.to('Myr') - get_time_offset(star_mode))**(2/5)
        if (st_radius < BOXSIZE_MIN):
            radius = radius_min
        else :
            radius = st_radius
    sphere = ds.sphere(center, radius)
    node.sphere = sphere


def retrieve_positions(node, position_file):
    df = h5py.File(position_file, 'r')
    ds = node.ds
    dumps = [dump.astype(str).split('/')[0] for dump in list(df['dumps'])]
    index = np.argwhere(np.array(dumps) == str(ds)).item()
    star_position = ds.arr(df['star_positions'][index], 'code_length').to('pc')
    halo_position = ds.arr(df['halo_positions'][index], 'code_length').to('pc')
    distance = unyt_quantity(np.linalg.norm(star_position - halo_position), 'pc')
    node.star_position = star_position
    node.halo_position = halo_position
    node.distance = distance


def align_sphere(node):
    ds = node.ds
    sp0 = node.sphere
    bulk_vel = sp0.quantities.bulk_velocity()
    sphere = ds.sphere(sp0.center, sp0.radius)
    sphere.set_field_parameter("bulk_velocity", bulk_vel)
    virial_radius = sphere.radius.in_units('pc').value.item()
    L = [sphere.quantities.total_quantity(('gas', f'angular_momentum_{ax}')) for ax in 'xyz']
    L_norm = [float((ax / np.linalg.norm(L)).value) for ax in L]
    sphere.set_field_parameter("normal", L_norm)
    node.sphere = sphere


def reunit(ds, val, units):
    if isinstance(val, yt.YTQuantity):
        func = ds.quan
    else:
        func = ds.arr
    return func(val.to(units).d, units)


def uperiodic_distance(x1, x2, domain=None):
    units = x1.units
    x2.convert_to_units(units)
    if domain is None:
        dom = None
    else:
        dom = domain.to(units).v

    d = periodic_distance(x1.v, x2.v, dom)
    return d * x1.uq


def modify_grackle(node, **grackle_params):

    ds = node.ds
    cooling_file = os.path.basename(ds.parameters["CloudyCoolingGridFile"])
    cooling_data_dir = "/home/brs/cloudy/cloudy_cooling_data"
    grackle_data_file = bytes(os.path.join(cooling_data_dir, cooling_file), 'utf-8')
    
    grackle_pars = {'grackle_data_file': grackle_data_file, 'cmb_temperature_floor': 0}
    for key, value in grackle_params.items():
        grackle_pars[key] = value
    pygrackle.add_grackle_fields(ds, parameters=grackle_pars)

    #print (sp['gas', 'grackle_cooling_time'])
    #print (sp['gas', 'grackle_gamma'])
    #print (sp['gas', 'grackle_mean_molecular_weight'])
    #print (sp['gas', 'grackle_pressure'])
    #print (sp['gas', 'grackle_temperature'])
    #print (sp['gas', 'grackle_dust_temperature'])

