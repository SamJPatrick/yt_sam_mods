import os
import sys
import yt
import numpy as np

from unyt import kb, mh, G
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib.pyplot as plt


MU = 1.6
ALPHA = 2
GAMMA = 5/3

OUTDIR = '.'
PRE_DIR = "Britton_sim_data/Clump_data/Profiles/"


def calculate_cool_rate(n_dim, T_dim):

    n = n_dim.to('cm**(-3)').value
    T = T_dim.to('K').value
    #A = - 178.4239 - 68.42243 * np.log(T) + 43.20243 * np.log(T)**2 - 4.633167 * np.log(T)**3 + 69.70086 * np.log(1 + 4.087038*10**4 / T)
    #B = 1.288953*10**2 - 53.91334 * np.log(T) + 5.315517 * np.log(T)**2 - 19.73427 * np.log(1 + 1.678095*10**4 / T)
    #C_ln = 14.82123 - 4.890915 * np.log(T) + 0.4749030 * np.log(T)**2 - 0.01338283 / T
    #C = np.exp(C_ln)
    #D = 0.8227443 + 0.5864073 * np.exp(-T / 1850) - 2.056313 * np.exp(-T /440)
    #eta_ln = A - \
    #    ((A - B) / (1 + n / C)**D) -\
    #    (2.370570*10**4 / T) + \
    #    ((-2.370570*10**4 / T) - (2.578611*10**4 / T)) * (1 + (n / (C * np.exp(-1.164408)))**D)**(-1)
    eta_log = -103.0 + 97.59 * np.log10(T) - 48.05 * np.log10(T)**2 + 10.80 * np.log10(T)**3 - 0.9032 * np.log10(T)**4
    eta = unyt_quantity(10**(eta_log), 'erg*cm**6/s')
    cool_rate = eta * n_dim**2
    return cool_rate

    
try:
    star_type = sys.argv[1]
except IndexError:
    star_type = 'pisn'
if (star_type == 'pisn' or star_type == 'PISN'):
    ds = yt.load(os.path.join(PRE_DIR, "DD0295_profile_number_density_temperature_radius.h5"))
elif (star_type == 'ccsn' or star_type == 'CCSN'):
    ds = yt.load(os.path.join(PRE_DIR, "DD0179_profile_number_density_temperature_radius.h5"))
elif (star_type == 'hn' or star_type == 'HN'):
    ds = yt.load(os.path.join(PRE_DIR, "DD0232_profile_number_density_temperature_radius.h5"))
else :
    print("Error, star type must be PISN, HN or CCSN")
    quit()
    

number_densities = ds.data[('data', 'x')].to('cm**(-3)')
LEN_N = len(number_densities)
temperatures = ds.data[('data', 'y')]
LEN_T = len(temperatures)
radii = ds.data[('data', 'z')]
LEN_R = len(radii)


argmin = 0
argmax = LEN_N
for i in range (LEN_N):
    if (1 in ds.data[('data', 'used')][i]):
        argmin = i
        break
for i in range (argmin, LEN_N):
    if (1 in ds.data[('data', 'used')][i]):
        continue
    else :
        argmax = i
        break
nmin = number_densities[argmin]
nmax = number_densities[argmax]


#total_mass = np.sum([np.sum([ds.data[('data', 'cell_mass')][n][:]]) for n in range (LEN_N)])
m_enc = unyt_array([np.nansum(
    [np.nansum(ds.data[('data', 'cell_mass')][i].to('Msun')) for i in range (n, LEN_N)])
                    for n in range (0, LEN_N)], 'Msun')
mass = unyt_array([np.nansum(
    ds.data[('data', 'cell_mass')][i]).to('Msun') for i in range (LEN_N)], 'Msun')
vol_enc = unyt_array([np.nansum(
    [(np.nansum(ds.data[('data', 'cell_mass')][i]) / number_densities[i]).to('g*cm**3') for i in range (n, LEN_N)])
                      for n in range (0, LEN_N)], 'g*cm**3')
n_ave = transpose_unyt([m / V for m, V in zip(m_enc, vol_enc)]).to('cm**(-3)')
r_enc = unyt_array([np.nansum(
    [(np.nansum(transpose_unyt([ds.data[('data', 'cell_mass')][i][k][j] for k in range (LEN_T)])) * radii[j]).to('g*m') for j in range (LEN_R)])
                    for i in range (LEN_N)], 'g*m')
r_ave = transpose_unyt([r / m for r, m in zip(r_enc, mass)]).to('pc')
#r0 = r_ave[argmin]

T_enc = unyt_array([np.nansum(
    [np.nansum(
        [(np.nansum(ds.data[('data', 'cell_mass')][i][j]) * temperatures[j]).to('g*K') for j in range (LEN_T)])
     for i in range (n, LEN_N)])
                    for n in range (0, LEN_N)], 'g*K')
T_ave = transpose_unyt([T / m for T, m in zip(T_enc, m_enc)]).to('K')


energy_enc = ((kb * (3.0/2.0) / (MU * mh)) * T_enc).to('erg')
#tau_cool = transpose_unyt([E / calculate_cool_rate(n, T) for E, n, T in zip(energy_enc, number_densities, T_ave)]).to('Myr')
tau_ff = np.sqrt(((3*np.pi)/(32*G))* (1/(MU * mh * n_ave))).to('Myr')
#tau_cs_theo = (r0 * (nmax / number_densities)**(1/ALPHA) * np.sqrt((MU * mh) / (GAMMA * kb * T_ave))).to('Myr')
tau_cs_exp = (r_ave * np.sqrt((MU * mh) / (GAMMA * kb * T_ave))).to('Myr')

index_co = np.nanargmin(np.abs(tau_ff - tau_cs_exp))
mass_be = (GAMMA / (MU * mh))**2 * (kb * T_ave[index_co])**(3/2) / np.sqrt(4 * np.pi * G**3 * number_densities[index_co])
mass_co = m_enc[index_co]
be_ratio = mass_co / mass_be

fpath = os.path.join(os.getcwd(), OUTDIR, f"{star_type}_clump_density_profiles.h5")
number_densities.in_units('cm**(-3)').write_hdf5(fpath, group_name='density')
r_ave.in_units('pc').write_hdf5(fpath, group_name='radius')
mass.in_units('Msun').write_hdf5(fpath, group_name='mass')
m_enc.in_units('Msun').write_hdf5(fpath, group_name='mass_enc')
T_ave.in_units('K').write_hdf5(fpath, group_name='temperature')
tau_ff.in_units('Myr').write_hdf5(fpath, group_name='free-fall')
tau_cs_exp.in_units('Myr').write_hdf5(fpath, group_name='sound-crossing')
energy_enc.in_units('erg').write_hdf5(fpath, group_name='energy_enc')
mass_be.in_units('Msun').write_hdf5(fpath, group_name='mass_be')
mass_co.in_units('Msun').write_hdf5(fpath, group_name='mass_co')
be_ratio.write_hdf5(fpath, group_name='be_ratio')
