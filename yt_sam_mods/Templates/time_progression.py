import numpy as np
import os
import yt
# yt.enable_parallelism()
import ytree
import matplotlib.pyplot as plt

from ytree.analysis import AnalysisPipeline

from yt.extensions.p2p.profiles import my_profile
from yt.extensions.p2p.tree_analysis_operations import \
    yt_dataset, \
    delattrs, \
    garbage_collect, \
    node_sphere

def file_does_not_exist(node, output_dir="."):
    # If the profile file has already been made, return False
    # so we skip it.
    fn = os.path.join(output_dir, f"node_{node.uid}_profiles.h5")
    return not os.path.exists(fn)

def metal_profile(node, output_dir="."):
    yt.mylog.info(f"Profiling {node}.")
    sphere = node.sphere
    profile = my_profile(sphere, [("index", "radius")], [("gas", "metallicity3")],
                         weight_field=("gas", "cell_mass"),
                         bin_density=10, units={("index", "radius"): "pc",
                                                ("gas", "metallicity3"): "Zsun"})
    fn = os.path.join(output_dir, f"node_{node.uid}_profiles.h5")
    profile.save_as_dataset(filename=fn)

def max_metallicity_in_sphere(node):
    yt.mylog.info(f"Checking metallicity of {node}.")
    sphere = node.sphere
    min_Z = sphere.ds.quan(1e-6, "Zsun")
    max_Z = sphere["gas", "metallicity3"].max()
    yt.mylog.info(f"{node} has maximum metallicity {max_Z}.")
    return max_Z >= min_Z

if __name__ == "__main__":
    output_dir = "."
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue"

    a = ytree.load(os.path.join(data_dir, "merger_trees/p2p_nd_fs/p2p_nd_fs.h5"))
    a.add_vector_field("icom_gas_position")
    my_uids = [43576459]
    my_trees = list(a[:])
    my_tree = a[0]
    node = my_tree.get_node("forest", 166118)

    # To get all the nodes:
    my_nodes = []
    while node is not None:
        my_nodes.append(node)
        node = node.descendent
        masses = a.arr([node["mass"] for node in my_nodes])
        
    metal_dir =	f"time_profiles/tree_{my_tree.uid}"
    ap = AnalysisPipeline(output_dir=output_dir)
    ap.add_operation(file_does_not_exist, output_dir=metal_dir)
    ap.add_operation(yt_dataset, data_dir)
    ap.add_operation(node_sphere, center_field="icom_gas_position")
    ap.add_operation(max_metallicity_in_sphere)
    ap.add_operation(metal_profile, output_dir=metal_dir)
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

   # prog =  list(my_tree["prog"])
    for node in my_nodes:
        print(metal_dir)
        ap.process_target(node)

            
       


