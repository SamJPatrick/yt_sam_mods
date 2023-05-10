import ytree

HALO_FILE = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"
POSITION_FILE = "halo_positions_test.txt"

a = ytree.load(HALO_FILE)
positions = list(zip(*[a[0]['prog', f'icom_gas_position_{ax}'].value.tolist() for ax in 'xyz']))
positions = [list(position) for position in positions]

with open(POSITION_FILE, 'w') as position_file:
    for position in positions:
        position_file.write(f"{position}\n")
