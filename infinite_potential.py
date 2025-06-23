#!/usr/bin/python3

# Schrödinger equation

from scipy.linalg import eigh_tridiagonal

import numpy as np

m = 9.10938356e-31  # electron mass in kg

L = 1e-10  # 1 angstrom in meters

N = 1000
dx = L / N  # spacial step

h = 6.62607015e-34  # Planck's constant in J.s
h_bar = h / (2 * np.pi)  # Reduced Planck's constant in J.s

main_diag = np.zeros(N)
sub_diag = -1.0 * np.ones(N - 1)


for k in range(N):
    x = k * dx
    vk = 0  # null potential
    bk = 2

    main_diag[k] = bk


eigh_values, eigh_vectors = eigh_tridiagonal(main_diag, sub_diag)

# converting eigenvalues to energy
energies = eigh_values * h_bar**2 / (2 * m * dx**2)  # converting to energy in Joules

# Converting energy values from Joules to meV
energies = energies / 1.602176634e-16  # J to meV

for i in range(19):
    print(f"State {i+1}: Energy = {energies[i]} meV")