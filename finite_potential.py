#!/usr/bin/python3

# Schrödinger equation

from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt


import numpy as np

m = 9.10938356e-31  # electron mass in kg

L_well = 1e-10  # 1 angstrom in meters
L_domain = 2*L_well  # 1 angstrom in meters

N = 1000
dx = L_domain / N  # spacial step

h = 6.62607015e-34  # Planck's constant in J.s
h_bar = h / (2 * np.pi)  # Reduced Planck's constant in J.s

const = h_bar**2/(2*m*dx**2)  # constant
main_diag = np.zeros(N)

v0 = 1.0e-19 # potential energy in J
#v0 = 0


l1_pot = -L_well/2  # start of the first well
l2_pot = L_well/2  # end of the first well

x_values = np.linspace(-L_domain/2, L_domain/2, N)
v = np.zeros(N)  # potential energy array


for i in range(N):

    if x_values[i] < l1_pot or x_values[i] > l2_pot:
        v[i] = v0
    else:
        v[i] = 0
        

main_diag = 2*const * np.ones(N) + v
sub_diag = -const * np.ones(N - 1)

eigh_values, eigh_vectors = eigh_tridiagonal(main_diag, sub_diag)

# converting eigenvalues to energy
energies = eigh_values * h_bar**2 / (2 * m * dx**2)  # converting to energy in Joules

# Converting energy values from Joules to eV
energies = energies / 1.602176634e-19  # J to eV


for i in range(19):
    print(f"State {i+1}: Energy = {energies[i]} eV")