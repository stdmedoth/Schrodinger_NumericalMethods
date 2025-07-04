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

const = - h_bar**2/(2*m*dx**2)  # constant

l1_pot = -L_well/2  # start of the first well
l2_pot = L_well/2  # end of the first well

x_values = np.linspace(-L_domain/2, L_domain/2, N)
v = np.zeros(N)  # potential energy array

v0 = 1.0e-10 # potential energy in J

for i in range(N):

    if x_values[i] < l1_pot or x_values[i] > l2_pot:
        v[i] = v0
    else:
        v[i] = 0
        

main_diag = 2*const * np.ones(N) + v
sub_diag = -const * np.ones(N - 1)

eigh_values, eigh_vectors = eigh_tridiagonal(main_diag, sub_diag)


#plt.plot(x_values, eigh_vectors[:, 0]**2, label='First ground state', color='blue')
#plt.plot(x_values, eigh_vectors[:, 1]**2, label='Second ground state', color='red')
plt.plot(x_values, eigh_vectors[:, 2]**2, label='Third ground state', color='green')


plt.axvline(l1_pot, color='red', lw=0.5, ls='--')  # horizontal line at y=0
plt.axvline(l2_pot, color='red', lw=0.5, ls='--')  # horizontal line at y=0

if np.max(v) > 0:
    # A escala é feita em relação ao valor máximo da primeira função de onda.
    # É mais comum usar o estado fundamental (índice 0) para essa normalização visual.
    scaling_factor = np.max(eigh_vectors[:, 0]**2) / np.max(v)
    plt.plot(x_values, v * scaling_factor, color='black', label='Scaled Potencial')
else:
    # Se o potencial for todo zero (poço infinito), plote-o sem escala, ou simplesmente não plote se não for relevante
    # Para o caso de v0=0, o potencial é uma linha em y=0. Podemos plotar isso diretamente.
    plt.plot(x_values, v, color='black', label='Potencial')


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

