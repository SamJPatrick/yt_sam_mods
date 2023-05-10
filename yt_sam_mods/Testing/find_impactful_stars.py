import os
import yaml
import numpy as np
import yt
import ytree
from unyt import Msun, Zsun, Myr
from unyt import unyt_array, unyt_quantity
from yt.extensions.sam_mods.unyt_funcs import *
from yt.extensions.sam_mods.misc import transpose_unyt




SIM_DIR = '/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue'
OUT_DIR = '.'




if __name__ == "__main__":

    with open(os.path.join(SIM_DIR, 'simulation.yaml')) as stream:
        snap_list = yaml.load(stream, Loader=yaml.FullLoader)
        snap_dict = dict([(entry['Snap_idx'], entry['filename'].split('/')[-1]) for entry in snap_list])

    with open(os.path.join(SIM_DIR, 'star_hosts.yaml')) as stream:
        df_halo = yaml.load(stream, Loader=yaml.FullLoader)

    sim = yt.load(os.path.join(SIM_DIR, 'simulation.h5'))
    efns = sim.data['filename'].astype(str)
    etimes = sim.data['time'].to('Myr')

    pop3_dir = os.path.join(SIM_DIR, 'pop3')
    info_dict = {}
    pop3_past = {}
    halos_past = {}
    for idx, dump_num in sorted(snap_dict.items()):
        
        time = etimes[np.argwhere(efns == f'{dump_num}/{dump_num}')].item() * Myr
        df_pop3 = yt.load(os.path.join(pop3_dir, f'{dump_num}.h5'))
        n_pop3 = len(df_pop3.data[('pop3', 'particle_index')])
        if (n_pop3 == 0):
            continue
        
        pop3_dict = {}
        halo_dict = {}
        for i in range (n_pop3):
            
            index = int(df_pop3.data[('pop3', 'particle_index')][i].item())
            position = [df_pop3.data[('pop3', f'particle_position_{ax}')][i].to('unitary') for ax in 'xyz']
            position = transpose_unyt(position)
            ct = df_pop3.data[('pop3', 'creation_time')][i].to('Myr')
            mass = df_pop3.data[('pop3', 'particle_mass')][i].to('Msun')
            if (mass < 1.0 * Msun):
                mass = pop3_past[index]['mass']
                if (pop3_past[index]['supernova_time'] == 0.0 * Myr):
                    st = time
                else :
                    st = pop3_past[index]['supernova_time']
            else :
                st = 0.0 * Myr
            pop3_past[index] = {'mass': mass, 'creation_time': ct, 'supernova_time': st, 'position': position}
            pop3_dict[index] = {'mass': unyt_quantity_to_yaml(mass), 'creation_time': unyt_quantity_to_yaml(ct),
                                'supernova_time': unyt_quantity_to_yaml(st), 'position': unyt_array_to_yaml(position)}

            if (index in halos_past.keys()):
                halo_dict[index] = {'mass': halos_past[index]['mass'], 'radius': halos_past[index]['radius'],
                                    'position': halos_past[index]['position']}
                continue
            halo_file = os.path.join(SIM_DIR, df_halo[index]['arbor'])
            tree_num = df_halo[index]['_arbor_index']
            node_num = df_halo[index]['tree_id']
            a = ytree.load(halo_file)
            node = a[tree_num].get_node('forest', node_num)
            mass = node['mass'].to('Msun')
            radius = node['virial_radius'].to('pc')
            position = node['position'].to('unitary')
            halo_info = {'mass': unyt_quantity_to_yaml(mass), 'radius': unyt_quantity_to_yaml(radius),
                         'position': unyt_array_to_yaml(position)}
            halo_dict[index] = halos_past[index] = halo_info
            
        info_dict[idx] = {'filename': dump_num, 'time': unyt_quantity_to_yaml(time), 'star_info': pop3_dict, 'halo_info': halo_dict}

    print("Finished building dictionary, outputing to yaml file 'star_halo_info.yaml'")

    with open(os.path.join(OUT_DIR, 'star_halo_info.yaml'), mode='wb') as stream:
        yaml.dump(info_dict, stream, encoding='utf-8')
