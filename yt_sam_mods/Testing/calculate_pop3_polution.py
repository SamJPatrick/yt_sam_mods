import os
import yaml
import numpy as np
import unyt
from unyt import Msun, Zsun, Myr, Mpc, erg
from unyt import unyt_array, unyt_quantity
from yt.extensions.sam_mods.unyt_funcs import *
from yt.extensions.sam_mods.misc import transpose_unyt




IN_DIR = '.'
OUT_DIR = '.'


BOX_SIZE = 0.5 * Mpc
RADIUS = unyt_quantity(1e3, 'pc')
POSITION = unyt_array([0.54673098, 0.52612445, 0.50558593], '')
#POSITION = unyt_array([16949.58529203, 18979.38793353, 17409.57691984], 'pc')



if __name__ == "__main__":

    with open(os.path.join(IN_DIR, 'star_halo_info.yaml')) as stream:
        info_dict = yaml.load(stream, Loader=yaml.FullLoader)

    supernova_dict = {}
    for idx_file, dp in info_dict.items():    
        assert sorted(dp['halo_info'].keys()) == sorted(dp['star_info'].keys()), \
            "Error, halo & star dictionaries must contain the same pop3 star keys!"
        
        for idx_star in dp['halo_info'].keys():
            # SHOULD POLUTE TIMES BE UPDATED????
            if (idx_star in supernova_dict.keys()):
                continue
            st = yaml_to_unyt_quantity(dp['star_info'][idx_star]['supernova_time']).to('Myr')
            if (st == 0.0 * Myr):
                continue
            
            radius = yaml_to_unyt_quantity(dp['halo_info'][idx_star]['radius'])
            position_star = yaml_to_unyt_array(dp['star_info'][idx_star]['position'])
            position_halo = yaml_to_unyt_array(dp['halo_info'][idx_star]['position'])
            offset = unyt_norm(position_star - position_halo) * BOX_SIZE
            if (offset  > radius):
                print(f"Error! Star is outside host halo for star {idx_star} in {dp['filename']}!!")
                print(f"Star-halo offset is {offset} whilst halo virial radius is {radius}") 
                continue
            
            mass_halo = yaml_to_unyt_quantity(dp['halo_info'][idx_star]['mass'])
            mass_star = yaml_to_unyt_quantity(dp['star_info'][idx_star]['mass'])
            rho = (3 * mass_halo) / (4 * np.pi * radius**3)
            distance = unyt_norm(position_star - POSITION) * BOX_SIZE
                
            if (mass_star < 140 * Msun):
                energy = unyt_quantity(1e51, 'erg')
                z_yeild = unyt_quantity(8, 'Msun')
            else :
                energy = unyt_quantity(1e53, 'erg')
                z_yeild = unyt_quantity(30, 'Msun')
            polute_time = (distance**(5.0/2.0) * (rho / energy)**(1.0/2.0)).to('Myr')
            z_reached = z_yeild * (RADIUS / distance)**2
            supernova_dict[idx_star] = {'polute_time': polute_time, 'metals': z_yeild, \
                                        'metals_frac': z_reached, 'supernova_time': st, \
                                        'time': yaml_to_unyt_quantity(dp['time'])}

    #supernova_dict_yaml = {}
    target_dict = {'n_polute': 0, 'tot_metals': 0, 'tot_metals_frac': 0}
    for idx_star, entry in supernova_dict.items():
        if (entry['time'] >= (entry['supernova_time'] + entry['polute_time'])):
            target_dict['n_polute'] += 1
            target_dict['tot_metals'] += entry['metals']
            target_dict['tot_metals_frac'] += entry['metals_frac']
        dictionary = supernova_dict[idx_star]
        for info, value in entry.items():
            dictionary[info] = unyt_quantity_to_yaml(value)
        #supernova_dict_yaml[idx_star] = dictionary

    with open(os.path.join(OUT_DIR, 'supernovae_info.yaml'), 'wb') as stream:
        yaml.dump(supernova_dict, stream, encoding='utf-8')

    with open(os.path.join(OUT_DIR, 'target_info.yaml'), 'wb') as stream:
        yaml.dump(target_dict, stream, encoding='utf-8')
