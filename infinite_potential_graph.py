#!/usr/bin/python3

# Schrödinger equation

from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt


import numpy as np

m = 9.10938356e-31  # electron mass in kg

L = 1e-10  # 1 angstrom in meters

N = 1000
dx = L / (N+1)  # spacial step

h = 6.62607015e-34  # Planck's constant in J.s
h_bar = h / (2 * np.pi)  # Reduced Planck's constant in J.s

const = - h_bar**2/(2*m*dx**2)  # constant
main_diag = np.zeros(N)
sub_diag = -1.0 * np.ones(N - 1) * const

for k in range(N):
    x = k * dx
    vk = 0 

    main_diag[k] = (-2*const) + vk


eigh_values, eigh_vectors = eigh_tridiagonal(main_diag, sub_diag)

x_values = np.linspace(0, L, N)

plt.plot(x_values, eigh_vectors[:, 0]**2, label='First ground state', color='blue')
#plt.plot(x_values, eigh_vectors[1]**2, label='Second ground state', color='red')
#plt.plot(x_values, eigh_vectors[2]**2, label='Third ground state', color='green')

plt.xlabel('Position (m)')
plt.ylabel('Probability Density')
plt.title('Ground State Probability Density in Infinite Potential Well')
plt.legend()
plt.grid()
plt.show()

# calculating the probability to find the electron in the middle of the well
probability_density = eigh_vectors[:, 0]**2
probability = 0

probability = np.sum(probability_density)
    
print(f"Probability to find the electron in the middle of the well = {probability:.4f}")

