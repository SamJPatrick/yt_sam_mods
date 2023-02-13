import re
import yt
import numpy as np
from unyt import unyt_quantity



SIMULATION = "simulation.h5"

T0_CCSNe = unyt_quantity(211.8104, 'Myr')
T0_PISNe = unyt_quantity(211.3103, 'Myr')

FIELD_DICTS = []
FIELD_DICTS.append({'name': 'density', 'units': 'g/cm**3', 'limits': (1e-25, 1e-18), 'log': True})
FIELD_DICTS.append({'name': 'temperature', 'units': 'K', 'limits': (1e1, 1e4), 'log': True})
FIELD_DICTS.append({'name': 'bonnor_ebert_ratio', 'units': None, 'limits': (1e-4, 1e2), 'log': True})
FIELD_DICTS.append({'name': 'velocity_magnitude', 'units': 'km/s', 'limits': (0, 1e6), 'log': False})
FIELD_DICTS.append({'name': 'sound_speed', 'units': 'km/s', 'limits': (0, 1e6), 'log': False})
FIELD_DICTS.append({'name': 'accretion_rate', 'units': 'Msun/yr', 'limits': (1e-5, 1e-2), 'log': True})
FIELD_DICTS.append({'name': 'accretion_rate_z', 'units': 'Msun/yr', 'limits': (1e-8, 1e-5), 'log': True})
FIELD_DICTS.append({'name': 'metallicity3', 'units': 'Zsun', 'limits': (1e-7, 1e0), 'log': True})


def get_field_dict(field):
    entry_num = np.argwhere(np.array([field['name'] for field in FIELD_DICTS]) == field)
    if (len(entry_num) != 1):
        print(f"Error, field {field} was not found within field list")
        quit()
    return FIELD_DICTS[entry_num.item()]


def get_time_offset(star_mode):
    time_offset = unyt_quantity(0.0, 'Myr')
    if (("PI" in star_mode) or ("pi" in star_mode)):
        time_offset = T0_PISNe
    elif (("CC" in star_mode) or ("cc" in star_mode)):
        time_offset = T0_CCSNe
    else:
        print("No offset provided, setting to zero...")
        pass
    return time_offset


def get_dump_num(filename):
    dump_num = re.search(r'^DD([0-9]{4})(_[-_a-zA-Z0-9]+\.h5)?$', filename).group(1)
    return dump_num


def get_time_z(filename, star_mode):
    time_offset = get_time_offset(star_mode)
    dump = f"DD{get_dump_num(filename)}"
    fname = '/'.join([dump, dump]).encode()
    sim = yt.load(SIMULATION)
    time = sim.data['time'][np.argwhere(sim.data['filename'] == fname).item()].in_units('Myr') - time_offset
    z = sim.data['redshift'][np.where(sim.data['filename'] == fname)][0].value
    return (time, z)


def get_title(filename, star_mode):
    time, z = get_time_z(filename, star_mode)
    dump = f"DD{get_dump_num(filename)}"
    #title = f"{dump}, z={z:.2f}, t = {time:.2f}"
    title = f"z={z:.2f}, t = {time:.2f}"
    return title
