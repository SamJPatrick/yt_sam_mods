import h5py
import glob
import os


DUMP_PATH = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/"
POSITIONS_FILE = "positions_info.h5"
HALO_FILE = "halo_positions_test.txt"
STAR_FILE = "pisn_positions_test.txt" 

df = h5py.File(POSITIONS_FILE, 'w')
with open(HALO_FILE, 'r') as halo_file:
    halo_positions = [[float(x) for x in arr.rstrip(']\n').lstrip('[').split(', ')] \
                      for arr in halo_file.readlines()]
with open(STAR_FILE, 'r') as star_file:
    star_positions = [[float(x) for x in arr.rstrip(']\n').lstrip('[').split(', ')] \
                     for arr in star_file.readlines()]
    
dumps = [os.path.basename(x) for x in glob.glob(os.path.join(DUMP_PATH, "DD*"))]
dumps = sorted(dumps, key= lambda x: int(x.lstrip('DD')))
dumps = [f'{dump}/{dump}'.encode("ascii", "ignore") for dump in dumps]

df.create_dataset("halo_positions", data= halo_positions)
df.create_dataset("star_positions", data= star_positions)
df.create_dataset("dumps", data= dumps)

