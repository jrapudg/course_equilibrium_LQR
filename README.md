# Course Equilibrium LQR

This project implements equilibrium state computation and LQR control for a wheeled robot model (Pineapple) using Julia. It uses the PinnZoo library for robot dynamics and kinematics, and implements optimization-based methods to find equilibrium states with contact constraints.

## Features

- Equilibrium state computation for wheeled robots with contact constraints
- KKT-based optimization using Newton's method
- Visualization using MeshCatMechanisms

## Installation

### Prerequisites

- Julia (version 1.6 or later)
- Git

### Setup Steps

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd course_equilibrium_LQR
   ```

2. **Initialize and update submodules:**
   ```bash
   git submodule update --init --recursive
   ```

3. **Switch to the pineapple branch in PinnZoo:**
   ```bash
   cd dependencies/PinnZoo
   git checkout pineapple
   ```

4. **Install PinnZoo**
   
    - Follow the README inside Pinnzoo folder
    - For Mac users, comment rigidbody project from `CMakeLists.txt` before installing PinnZoo

## Project Structure

- `LQR.jl` - Main script for equilibrium computation and LQR setup
- `dynamics_pineapple.jl` - Dynamics and contact constraint utilities for the Pineapple model
- `dependencies/PinnZoo/` - PinnZoo library submodule (must be on `pineapple` branch)

## Notes

- The project requires PinnZoo to be on the `pineapple` branch
- The model uses quaternion representation for orientation
- Contact constraints are enforced for wheel-ground contact
