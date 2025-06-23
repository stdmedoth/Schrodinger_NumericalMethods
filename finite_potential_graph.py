#!/usr/bin/python3

# Schrödinger equation

from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt


import numpy as np

m = 9.10938356e-31  # electron mass in kg

L = 1e-10  # 1 angstrom in meters

N = 2000
dx = L / N  # spacial step

h = 6.62607015e-34  # Planck's constant in J.s
h_bar = h / (2 * np.pi)  # Reduced Planck's constant in J.s

const = h_bar**2/(2*m*dx**2)  # constant
main_diag = np.zeros(N)
sub_diag = -const * np.ones(N - 1)

v0 = 5.0e-17 # potential energy in J
#v0 = 0


l1 = -L/2  # start of the first well
l2 = L/2  # end of the first well

xk = -L
for k in range(N):
    xk += k * dx
    
    vk = 0  # null potential outside the well
    if xk < l1 or xk > l2:
        vk = v0  # constant potential energy inside the well
    
    main_diag[k] = 2*const + vk


eigh_values, eigh_vectors = eigh_tridiagonal(main_diag, sub_diag)

x_values = np.linspace(-L, L, N)

plt.figure(figsize=(10, 6))

plt.axvline(l1, color='red', lw=0.5, ls='--')  # horizontal line at y=0
plt.axvline(l2, color='red', lw=0.5, ls='--')  # horizontal line at y=0


#plt.plot(x_values, eigh_vectors[:, 0]**2, label='First ground state', color='blue')
plt.plot(x_values, eigh_vectors[1]**2, label='Second ground state', color='green')
#plt.plot(x_values, eigh_vectors[2]**2, label='Third ground state', color='orange')

plt.xlabel('Position (m)')
plt.ylabel('Probability Density')
plt.title('Ground State Probability Density in finite Potential Well')
plt.legend()
plt.grid()
plt.show()

