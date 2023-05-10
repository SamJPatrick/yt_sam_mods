import yt
import os
import glob
import numpy as np
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt


STAR_DIR = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/pop3"
POSITION_FILE = "/disk01/spatrick/PISN_data/pisn_positions_test.txt"


datasets = sorted(glob.glob(os.path.join(STAR_DIR, "DD*.h5")))
ds = yt.load(datasets[0])
index_false = int(ds.data[('pop3', 'particle_index')])
ds = yt.load(datasets[1])
indices = [int(x) for x in ds.data[('pop3', 'particle_index')] if x != index_false]
assert len(indices) == 1, "Error, multiple stars created in 2nd dataset, correct index ambiguous"
index = indices[0]

positions = []
for dataset in datasets[1:]:
    ds = yt.load(dataset)
    indices = np.array([int(x) for x in ds.data[('pop3', 'particle_index')]])
    try :
        n = np.argwhere(indices == index).item()
    except IndexError:
        continue
    position = transpose_unyt([ds.data[('pop3', f'particle_position_{ax}')][n] for ax in 'xyz']).value.tolist()
    positions.append(position)

with open(POSITION_FILE, 'w') as position_file:
    for position in positions:
        position_file.write(f"{position}\n")
