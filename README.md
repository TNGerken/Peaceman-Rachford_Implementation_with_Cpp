# Peaceman-Rachford_Implementation_with_Cpp
Implementation of the Peaceman–Rachford scheme for solving transient heat diffusion PDEs. Includes detailed derivations, discretization steps, and C++ code for simulation and visualization.

# Peaceman–Rachford Heat Diffusion Solver

## Overview
This repository implements the Peaceman–Rachford method to solve transient heat diffusion problems in two dimensions. The project includes:

- Step-by-step derivation of the scheme
- Finite difference discretization
- Implementation in Python
- Example simulation with visual outputs

The Peaceman–Rachford scheme is an unconditionally stable method, making it well-suited for parabolic PDEs such as the heat equation, especially when high temporal accuracy and stability are desired.

---

## Methodology
The project solves the heat equation:

$$
\frac{\partial u}{\partial t} = \alpha \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)$$

using the Peaceman–Rachford ADI method, which splits each time step into two implicit half-steps:

1. Implicit in \(x\), explicit in \(y\)
2. Explicit in \(x\), implicit in \(y\)

This reduces the multidimensional implicit problem into simpler 1D tridiagonal systems, allowing efficient computation.

---

## Features
- Unconditionally stable time integration
- Modular Python code
- Easily adjustable parameters (grid size, time step, diffusivity)
- Visualization of temperature field evolution

---

## File Structure
Peaceman-Rachford_Implementation_with_Cpp/ 
│── data/ # Generated data for different timesteps
│── src/ # C++ main script
└── Documentation/ # Complete detailed report

---

## Features
- C++ implementation for computational efficiency.
- Handles large time steps with improved stability.
- Outputs numerical results to `data/` for analysis.
- Includes a detailed report explaining the method.

---


