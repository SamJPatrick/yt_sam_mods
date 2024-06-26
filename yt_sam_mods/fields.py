###Pop2Prime fields.

import numpy as np
from yt.utilities.exceptions import YTFieldNotFound
from unyt import kb, G, Msun, pc, mh
from unyt import unyt_quantity, unyt_array



def _vrot(field, data):
    v_rot = data[('gas', 'specific_angular_momentum_magnitude')] / data['radius']
    np.nan_to_num(v_rot, posinf=0.0, copy=False)
    return v_rot

def _overlap(field, data):
    bary_grads = list(zip(*[data[('gas', f'density_gradient_{ax}')] for ax in 'xyz']))
    dm_grads = list(zip(*[data[('gas', f'dark_matter_density_gradient_{ax}')] for ax in 'xyz']))
    grads_intcp = [np.dot(bary_grad, dm_grad) for bary_grad, dm_grad in zip(bary_grads, dm_grads)]
    align = [grad / np.abs(grad_bary_mag * grad_dm_mag) for grad, grad_bary_mag, grad_dm_mag in \
             zip(grads_intcp, data[('gas', 'density_gradient_magnitude')], data[('gas', 'dark_matter_density_gradient_magnitude')])]
    return align

#def _relative_radial_velocity(field, data):
#    x = data.get_field_parameter("reference_point")
#    locs = zip(*[data[('gas', '')]])


'''
def _vturb(field, data):
    data._debug()
    try :
        center = data.get_field_parameter('center')
        bulk_vel = np.linalg.norm(data.get_field_parameter('bulk_velocity').in_units('km/s'))
        bulk_vel = unyt_quantity(bulk_vel, 'km/s')
        #vturb = np.sqrt(data[('gas', 'velocity_magnitude')]**2 - data[('gas', 'vrot')]**2 - bulk_vel**2)
        vturb = np.sqrt(data[('gas', 'velocity_magnitude')]**2 - data[('gas', 'radial_velocity')]**2 - bulk_vel**2)
        np.nan_to_num(vturb, nan=0.0, copy=False)
        return vturb
    except Exception as ex:
        print(f"Error in vturb function code of type {type(ex).__name__} and arguments {ex.args}")
        quit()
    
def _cs_quad(field, data):
    data._debug()
    try :
        return np.sqrt(data[('gas', 'sound_speed')]**2 + data[('gas', 'vturb')]**2)
    except Exception as ex:
        print(f"Error in cs_quad function code of type {type(ex).__name__} and arguments {ex.args}")
        quit()
'''

        
def _metallicity3(field, data):
    return data[('enzo', 'SN_Colour')] / data[('gas', 'density')]

def _metallicity3_min7(field, data):
    field_data = data[('enzo', 'SN_Colour')] / data[('gas', 'density')]
    field_data.convert_to_units("")
    min_Z = unyt_quantity(1.e-7, "Zsun").in_units("")
    field_data[field_data < min_Z] = 0.5 * min_Z
    return field_data

def _metal3_mass(field, data):
    return data[('enzo', 'SN_Colour')] * data[('index', 'cell_volume')]

def _total_metal_density(field, data):
    field_data = np.zeros_like(data[('gas', 'density')])
    fields = [('enzo', 'Metal_Density'), ('enzo', 'SN_Colour')]
    for field in fields:
        if field in data.ds.field_list:
            field_data += data[field]
    return field_data

def _total_metallicity(field, data):
    return data[('gas', 'total_metal_density')] / data[('gas', 'density')]



# [erg] [m]^-3 ---> kinetic energy per unit volume
#def _kinetic_energy_density(field, data):
#    defined within yt

# [erg][kg]^-1 ---> thermal energy per unit masss
#def specific_thermal_energy(field, data):
#    defined within yt

# [erg][m]^-3 ---> total energy (thermal + kinetic) per unit volume
def _intrinsic_energy_density(field, data):
    return data[('gas', 'kinetic_energy_density')] + \
        data[('gas', 'specific_thermal_energy')] * data[('gas', 'density')]

# [erg] ---> total energy (thermal + kinetic)
def _intrinsic_energy(field, data):
    return data[('gas', 'intrinsic_energy_density')] * data[('gas', 'cell_volume')]

# [erg][kg]^-1 ---> kinetic energy per unit mass
def _specific_kinetic_energy(field, data):
    return data[('gas', 'kinetic_energy_density')] / data[('gas', 'density')]



def _pressure_accl(field, data):
    return - data[('gas', 'pressure_gradient_magnitude')] / data[('gas', 'density')]

def _inverse_jeans(field, data):
    return 1 / data[('gas', 'jeans_mass')]



def _l_phi(field, data):
    return data[('gas', 'velocity_spherical_phi')] * data[('gas', 'cell_mass')] * data['radius']

def _ang_ratio(field, data):
    return np.abs(data[('gas', 'velocity_spherical_phi')] / data[('gas', 'tangential_velocity')])

def _tangential_velocity_magnitude(field, data):
    return np.sqrt(data[('gas', 'velocity_spherical_theta')]**2 + data[('gas', 'velocity_spherical_phi')]**2)



def _mass_accretion_rate(field, data):
    return - 4 * np.pi * data[('index', 'radius')]**2 * data[('gas', 'density')] * data[('gas', 'velocity_spherical_radius')]

def _metal_mass_accretion_rate(field, data):
    return - 4 * np.pi * data[('index', 'radius')]**2 * data[('gas', 'density')] * data[('gas', 'metallicity3')] * data[('gas', 'velocity_spherical_radius')]



def _total_dynamical_time(field, data):
    return np.sqrt(3.0 * np.pi / (16.0 * data[('gas', 'matter_density')] * G))

def _sound_crossing_time(field, data):
    return data['radius'] / data[('gas', 'sound_speed')]

def _sound_crossing_time_mod(field, data):
    return data['radius'] / data[('gas', 'vsigma')]

def _cooling_ratio(field, data):
    return data['radius'] / (data[('gas', 'sound_speed')] * data[('gas', 'cooling_time')])

def _cooling_dynamical_ratio(field, data):
    return data[('gas', 'cooling_time')] / data[('gas', 'dynamical_time')]

def _cooling_rate(field, data):
    data._debug()
    return data[('gas', 'internal_energy')] * data[('gas', 'cell_mass')] / data[('gas', 'cooling_time')]

def _vortical_dynamical_ratio(field, data):
    return data[('gas', 'vortical_time')] / data[('gas', 'dynamical_time')]

def _vortical_time(field, data):
    return 1.0 / data[('gas', 'vorticity_magnitude')]

def _vortical_cooling_ratio(field, data):
    return data[('gas', 'vortical_time')] / data[('gas', 'cooling_time')]

def _cooling_intensity(field, data):
    n_density = data[('gas', 'density')] / (data[('gas', 'mean_molecular_weight')] * mh)
    return data[('gas', 'temperature')] * n_density**2



def _dark_matter_mass(field, data):
    return data[('gas', 'dark_matter_density')] * data[('index', 'cell_volume')]

def _HD_H2(field, data):
    return data[('gas', 'HD_density')] / data[('gas', 'H2_density')]






def add_p2p_field(ds, name, function=None, units='', sampling_type='cell'):
    try:
        ds.add_field(name, function=function,
                     units=units, sampling_type=sampling_type)
    except YTFieldNotFound:
        pass

def add_p2p_fields(ds):
    # use the value of solar metallicity in the dataset
    ds.unit_registry.modify('Zsun', ds.parameters['SolarMetalFractionByMass'])

    add_p2p_field(ds, ('gas', 'rotational_velocity'),
                  function=_vrot,
                  units='km/s', sampling_type='cell')
    #add_p2p_field(ds, ('gas', 'turbulent_velocity'),
    #              function=_vturb,
    #              units='km/s', sampling_type='cell')
    #add_p2p_field(ds, ('gas', 'cs_quad'),
    #              function=_cs_quad,
    #              units='km/s', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'metallicity3'),
                  function=_metallicity3,
                  units='Zsun', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'metallicity3_min7'),
                  function=_metallicity3_min7,
                  units='Zsun', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'metal3_mass'),
                  function=_metal3_mass,
                  units='Msun', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'total_metal_density'),
                  function=_total_metal_density,
                  units='g/cm**3', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'total_metallicity'),
                  function=_total_metal_density,
                  units='g/cm**3', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'intrinsic_energy'),
                  function=_intrinsic_energy,
                  units='erg', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'intrinsic_energy_density'),
                  function=_intrinsic_energy_density,
                  units='erg/m**3', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'specific_kinetic_energy'),
                  function=_specific_kinetic_energy,
                  units='erg/kg', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'pressure_accl'),
                  function=_pressure_accl,
                  units='km/s**2', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'inverse_jeans'),
                  function=_inverse_jeans,
                  units='1/Msun', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'l_phi'),
                  function=_l_phi,
                  units='pc*km/s', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'ang_ratio'),
                  function=_ang_ratio,
                  units='', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'tangential_velocity_magnitude'),
                  function=_tangential_velocity_magnitude,
                  units='km/s', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'accretion_rate'),
                  function=_mass_accretion_rate,
                  units='Msun/yr', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'accretion_rate_z'),
                  function=_metal_mass_accretion_rate,
                  units='Msun/yr', sampling_type='cell')

    add_p2p_field(ds, ('gas', 'total_dynamical_time'),
                  function=_total_dynamical_time,
                  units='yr', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'sound_time'),
                  function=_sound_crossing_time,
                  units='yr', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'sound_time_mod'),
                  function=_sound_crossing_time_mod,
                  units='yr', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'cooling_ratio'),
                  function=_cooling_ratio,
                  units='', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'cooling_dynamical_ratio'),
                  function=_cooling_dynamical_ratio,
                  units='', sampling_type='cell')
    #add_p2p_field(ds, ('gas', 'cooling_rate'),
    #              function=_cooling_rate,
    #              units='erg/s', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'vortical_dynamical_ratio'),
                  function=_vortical_dynamical_ratio,
                  units='', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'vortical_time'),
                  function=_vortical_time,
                  units='yr', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'vortical_cooling_ratio'),
                  function=_vortical_cooling_ratio,
                  units='', sampling_type='cell')
    add_p2p_field(ds, ('gas', 'cooling_intensity'),
                  function=_cooling_intensity,
                  units='K/cm**6', sampling_type='cell')


    
    add_p2p_field(ds, ('gas', 'dark_matter_mass'),
                  function=_dark_matter_mass,
                  units='Msun', sampling_type="cell")
    add_p2p_field(ds, ('gas', 'HD_H2_ratio'),
                  function=_HD_H2,
                  units='', sampling_type='cell')
    
    add_p2p_field(ds, ('gas', 'overlap'),
                  function=_overlap,
                  units='', sampling_type='cell')
