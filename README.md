# Quantum Box Numerical Schrodinger Solver

This repository contains Python code for numerically solving the Time-Independent Schrödinger Equation (TISE) for one-dimensional infinite and finite potential wells (also known as "particle in a box" problems).

The project implements a finite difference method to discretize the Schrödinger equation, transforming it into an eigenvalue problem that is solved using numerical diagonalization techniques from the `scipy.linalg` library.

## Features:
* Numerical solution for the infinite potential well.
* Numerical solution for the finite potential well.
* Calculation of energy eigenvalues (in Joules and electron-volts) and corresponding eigenfunctions.
* Utilizes `numpy` for efficient array operations and `scipy` for tridiagonal matrix diagonalization.

This work serves as a practical application of computational physics principles, aiming to provide a clear and reproducible approach to understanding quantum mechanical systems through numerical methods.

## Getting Started:
Clone the repository and run the Python scripts.
