import re
import yt
import numpy as np
from unyt import unyt_quantity



SIMULATION = "simulation.h5"

T0_CCSNe = unyt_quantity(211.31, 'Myr')
T0_HNe = unyt_quantity(211.31, 'Myr')
T0_PISNe = unyt_quantity(211.06, 'Myr')
#T0_PISNe = unyt_quantity(211.16, 'Myr')

TLIFE_CCSNe = unyt_quantity(3.86, 'Myr')
TLIFE_HNe = unyt_quantity(3.86, 'Myr')
TLIFE_PISNe = unyt_quantity(2.20, 'Myr')

E_CCSNe = unyt_quantity(1e51, 'erg')
E_HNe = unyt_quantity(3e52, 'erg')
E_PISNe = unyt_quantity(1e53, 'erg')

LIMS_MASS = (1e-8, 1e6)
LIMS_DENSITY = (1e-27, 1e-18)
LIMS_TEMPERATURE = (1e1, 1e8)
LIMS_PRESSURE = (1e-15, 1e-9)
LIMS_ENTROPY = (1e-2, 5e2)
LIMS_INTENSITY = (1e-6, 1e7)
LIMS_CTIME = (1e-4, 1e9)
LIMS_BE_RATIO = (1e-4, 1e2)
LIMS_VELOCITY = (-1e1, 3e1)
LIMS_CS = (0, 1e6)
LIMS_MDOT = (1e-5, 1e-2)
LIMS_MDOT_Z = (1e-8, 1e-5)
LIMS_Z = (1e-7, 1e0)
LIMS_H2 = (1e-4, 5e-3)
LIMS_EL = (1e-8, 1e0)

FIELD_DICTS = []
FIELD_DICTS.append({'name': 'cell_mass', 'units': 'Msun', 'limits': LIMS_MASS, 'colormap': 'viridis', 'log': True})
FIELD_DICTS.append({'name': 'density', 'units': 'g/cm**3', 'limits': LIMS_DENSITY, 'colormap': 'viridis', 'log': True})
FIELD_DICTS.append({'name': 'dark_matter_density', 'units': 'g/cm**3', 'limits': LIMS_DENSITY, 'colormap': 'viridis', 'log': True})
FIELD_DICTS.append({'name': 'temperature', 'units': 'K', 'limits': LIMS_TEMPERATURE, 'colormap': 'plasma', 'log': True})
FIELD_DICTS.append({'name': 'pressure', 'units': 'Pa', 'limits': LIMS_PRESSURE, 'colormap': 'spring', 'log': True})
FIELD_DICTS.append({'name': 'entropy', 'units': 'erg*cm**2', 'limits': LIMS_ENTROPY, 'colormap': 'spring', 'log': True})
FIELD_DICTS.append({'name': 'cooling_intensity', 'units': 'K/cm**6', 'limits': LIMS_INTENSITY, 'colormap': 'plasma', 'log': True})
FIELD_DICTS.append({'name': 'cooling_time', 'units': 'Myr', 'limits': LIMS_CTIME, 'colormap': 'plasma', 'log': True})
FIELD_DICTS.append({'name': 'bonnor_ebert_ratio', 'units': '', 'limits': LIMS_BE_RATIO, 'colormap': 'cividis', 'log': True})
FIELD_DICTS.append({'name': 'velocity', 'units': 'cm/s', 'limits': LIMS_VELOCITY, 'colormap': 'spring', 'log': False}) # <-----  Change back to 'km/s' once done
FIELD_DICTS.append({'name': 'sound_speed', 'units': 'km/s', 'limits': LIMS_CS, 'colormap': 'spring', 'log': False})
FIELD_DICTS.append({'name': 'accretion_rate', 'units': 'Msun/yr', 'limits': LIMS_MDOT, 'colormap': 'magma', 'log': True})
FIELD_DICTS.append({'name': 'accretion_rate_z', 'units': 'Msun/yr', 'limits': LIMS_MDOT_Z, 'colormap': 'magma', 'log': True})
FIELD_DICTS.append({'name': 'metallicity3', 'units': 'Zsun', 'limits': LIMS_Z, 'colormap': 'cool', 'log': True})
FIELD_DICTS.append({'name': 'H2_p0_fraction', 'units': '', 'limits': LIMS_H2, 'colormap': 'cool', 'log': False})
FIELD_DICTS.append({'name': 'El_fraction', 'units': '', 'limits': LIMS_EL, 'colormap': 'cool', 'log': True})


def get_field_dict(field):
    if (isinstance(field, str)):
        my_field = field
    elif (isinstance(field, tuple)):
        my_field = field[1]
    else :
        raise TypeError(f"Error, expecting either tuple or string for field, instead got {field}")
    
    if ("velocity" in my_field):
        entry_num = np.argwhere(np.array([field['name'] for field in FIELD_DICTS]) == "velocity")
    else :
        entry_num = np.argwhere(np.array([field['name'] for field in FIELD_DICTS]) == my_field)
    if (len(entry_num) != 1):
        print(f"Error, field {my_field} was not found within field list")
        return None
    
    return FIELD_DICTS[entry_num.item()]



def get_time_offset(star_mode):
    time_offset = unyt_quantity(0.0, 'Myr')
    if (("PI" in star_mode) or ("pi" in star_mode)):
        time_offset = T0_PISNe
    elif (("CC" in star_mode) or ("cc" in star_mode)):
        time_offset = T0_CCSNe
    elif (("HN" in star_mode) or ("hn" in star_mode)):
        time_offset = T0_HNe
    else:
        print("No offset provided, setting to zero...")
        pass
    return time_offset


def get_lifetime_offset(star_mode):
    time_offset = unyt_quantity(0.0, 'Myr')
    if (("PI" in star_mode) or ("pi" in star_mode)):
        time_offset = TLIFE_PISNe
    elif (("CC" in star_mode) or ("cc" in star_mode)):
        time_offset = TLIFE_CCSNe
    elif (("HN" in star_mode) or ("hn" in star_mode)):
        time_offset = TLIFE_HNe
    else:
        print("No offset provided, setting to zero...")
        pass
    return time_offset


def get_sn_energy(star_mode):
    if (("PI" in star_mode) or ("pi" in star_mode)):
        sn_energy = E_PISNe
    elif (("CC" in star_mode) or ("cc" in star_mode)):
        sn_energy = E_CCSNe
    elif (("HN" in star_mode) or ("hn" in star_mode)):
        sn_energy = E_HNe
    else:
        sn_energy = E_CCSNe
        print("None provided, setting to values for core collapse...")
        pass
    #print(f"Energy set to: {sn_energy}") 
    return sn_energy


def get_dump_num(filename):
    dump_num = re.search(r'^DD([0-9]{4})(_[-_a-zA-Z0-9]+\.h5)?$', filename).group(1)
    return dump_num


def get_time_z(filename, star_mode, sim_file=None):
    time_offset = get_time_offset(star_mode)
    dump = f"DD{get_dump_num(filename)}"
    fname = '/'.join([dump, dump]).encode()
    if (sim_file == None):
        sim = yt.load(SIMULATION)
    else :
        sim = yt.load(sim_file)
    time = sim.data['time'][np.argwhere(sim.data['filename'] == fname).item()].in_units('Myr') - time_offset
    z = sim.data['redshift'][np.where(sim.data['filename'] == fname)][0].value
    return (time, z)


def get_title(filename, star_mode):
    time, z = get_time_z(filename, star_mode)
    dump = f"DD{get_dump_num(filename)}"
    #title = f"{dump}, z={z:.2f}, t = {time:.2f}"
    title = f"z={z:.2f}, t = {time:.2f}"
    return title
