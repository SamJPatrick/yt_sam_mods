import yt
import pygrackle
import os
import numpy as np
import statistics
from unyt import unyt_array, unyt_quantity
from unyt import pc, Myr
from yt.extensions.astro_analysis.halo_analysis.halo_catalog.halo_callbacks import periodic_distance
from yt.utilities.exceptions import YTSphereTooSmall




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
    sphere = ds.sphere(center, radius)
    node.sphere = sphere


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


def transpose_unyt(quant_arr):
    values, units = zip(*[(quantity.value.item(), quantity.units) for quantity in quant_arr])
    unit = statistics.mode(units)
    arr = unyt_array([quantity.in_units(unit).value.item() for quantity in quant_arr], unit)
    return arr


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

