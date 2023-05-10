import h5py
import os
import numpy as np
from unyt import unyt_array, unyt_quantity, pc


DIR = "/disk01/spatrick/PISN_data/Sobolev/Ray_profiles"
FILE = "DD0160_sobolev.h5"
NUM_RAYS = 10
# both values in pc
BUMP_GUESS = 270 * pc
WINDOW = 60 * pc
MAG = 1.5


print("Ray stats")
print(f"Local maximum: ~ {BUMP_GUESS}")
print(f"Window width: {WINDOW}")
print(f"Magnitude drop: {MAG}")
print("\n")


df = h5py.File(os.path.join(DIR, FILE), 'r')
lengths = []
for n in range (NUM_RAYS - 1):
    densities = unyt_array(df[f'density_{n+1}/array_data'], 'g/cm**3')
    radii = unyt_array(df[f'distances_{n+1}/array_data'], 'pc')

    close = (np.abs(radii - BUMP_GUESS) < WINDOW)
    high = (densities * 10**MAG > np.max(densities[close]))
    bump_mask = (close & high)
    index_max = np.argwhere(densities == np.max(densities[close])).item()

    print(f"Stats for ray # {n+1}")
    print("Radii around local maximum")
    print(radii[bump_mask])
    print("Radius of local maximum")
    print(radii[index_max])
    print("Density of local maximum")
    print(np.max(densities[close]))

    if (len(radii[bump_mask]) <= 2):
        print("\n")
        continue
    if (radii[index_max] == radii[bump_mask][0]):
        radius_max = radii[bump_mask][1]
    elif (radii[index_max] == radii[bump_mask][-1]):
        radius_max = radii[bump_mask][-2]
    else :
        radius_max = radii[index_max]

    grad_1 = (np.max(densities[close]) - densities[bump_mask][0]) / (radius_max - radii[bump_mask][0])
    grad_2 = (densities[bump_mask][-1] - np.max(densities[close])) / (radii[bump_mask][-1] - radius_max)
    grad_grad = (grad_2 - grad_1) / (radii[bump_mask][-1] - radii[bump_mask][0])
    l_sob = np.sqrt(np.max(densities[close]) / np.abs(grad_grad))
    print(l_sob)
    print("\n")
    lengths.append(l_sob)

l_mean = np.mean(lengths)
print("\n")
print(f"Mean length scale: {l_mean:.2f}")
