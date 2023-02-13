import numpy as np
import os
import time
import yaml
import yt
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset

from ytree.analysis import AnalysisPipeline
from yt.extensions.p2p.misc import sphere_icom, reunit
from yt.extensions.p2p.profiles import my_profile
from yt.extensions.p2p.tree_analysis_operations import *
from yt.extensions.p2p import add_p2p_fields

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc




def get_dataset_node(ds, tree_path):

    a = ytree.load(tree_path)
    time = ds.current_time.in_units('Myr')
    node_times = [node['time'] for node in list(a[0]['prog'])]
    ifn = np.abs(node_times - time).argmin()
    node = list(a[0]['prog'])[ifn]
    return node





if __name__ == "__main__":

    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"
    dump_num = "DD0288/DD0288"

    dump_path = os.path.join(data_dir, dump_num)
    ds = yt.load(dump_path)
    node = get_dataset_node(ds, tree_path)
    sphere = get_node_sphere(node, ds= ds, center_field= 'position')

    my_profile(sphere, [('index', 'radius'), ('gas', 'vortical_time')], ('gas', 'cell_mass'), \
               weight_field= None, bin_density= 20, logs= {('index', 'radius'): True, ('gas', 'vortical_time'): True}, \
               units= {('index', 'radius'): 'pc', ('gas', 'vortical_time'): 'yr'})
    my_profile(sphere, [('index', 'radius')], ('gas', 'vortical_time'), weight_field= ('gas', 'cell_mass'), \
               bin_density= 20, logs= {('index', 'radius'): True, ('gas', 'vortical_time'): True}, \
               units= {('index', 'radius'): 'pc', ('gas', 'vortical_time'): 'yr'})
