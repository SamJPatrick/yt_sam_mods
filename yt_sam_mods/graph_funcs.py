import re
import yt
import numpy as np



sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
#sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"


def get_dump_num(filename):

    dump_num = re.search(r'^DD([0-9]{4})_[-_a-zA-Z]+.h5$', filename).group(1)
    return dump_num


def get_time_z(filename):
    
    dump = f"DD{get_dump_num(filename)}"
    fname = '/'.join([dump, dump]).encode()
    sim = yt.load(sim_path)
    time = sim.data['time'][np.argwhere(sim.data['filename'] == fname).item()].in_units('Myr')
    z = sim.data['redshift'][np.where(sim.data['filename'] == fname)][0].value
    return (time, z)


def get_title(filename):

    time, z = get_time_z(filename)
    dump = f"DD{get_dump_num(filename)}"
    title = f"{dump}, z={z:.2f}, t = {time:.2f}"
    return title
