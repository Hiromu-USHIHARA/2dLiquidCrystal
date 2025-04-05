
# 2D Liquid Crystal Simulation

This repository contains a Monte Carlo simulation of a 2D nematic liquid crystal system using the Metropolis algorithm. Molecules are modeled as rigid rods in 2D space, each with a position and an orientation. The simulation visualizes the evolution of the system under thermal fluctuations and optionally an external field.

<p align="center">
  <img src="LC.num256.time128.gif" alt="Liquid Crystal Simulation GIF" width="500">
</p>

---

## 📌 Features
- **2D liquid crystal simulation** using Metropolis Monte Carlo
- **Rods with continuous orientation angles**
- Interaction based on spatial and angular overlap
- Visualization via animated GIF and order parameter plot

---

## 📁 File Structure

| File                          | Description                                  |
|-------------------------------|----------------------------------------------|
| `nematic_monte_carlo.jl`     | Main simulation script                       |
| `LC.num256.time128.gif`      | Animated evolution of liquid crystal system  |
| `S.num256.time128.png`       | Order parameter $S(t)$ vs time plot          |
| `LC.data.num256.time128.txt` | Raw simulation data (positions and angles)   |

---

## ▶️ Usage

### 1. Prerequisites
Install Julia (≥1.6 recommended) and the required packages:

```julia
using Pkg
Pkg.add("Plots")
```

### 2. Run the Simulation
```bash
julia nematic_monte_carlo.jl
```

This will generate:
- `LC.num256.time128.gif`: animation of rod orientations
- `S.num256.time128.png`: order parameter over time
- `LC.data.num256.time128.txt`: raw system data

---

## 📊 Order Parameter
The simulation tracks the **nematic order parameter** $S$ over time:

$$
S = \sqrt{\langle \cos(2	heta) 
angle^2 + \langle \sin(2	heta) 
angle^2}
$$

This measures the degree of molecular alignment in the system.

<p align="center">
  <img src="S.num256.time128.png" alt="Order Parameter Plot" width="400">
</p>

---

## ⚙️ Adjustable Parameters
In `nematic_monte_carlo.jl`, modify the constants below to change simulation behavior:

```julia
const num_particles = 2^8        # Number of rod-like molecules
const num_steps = 2^7            # Simulation steps
const system_size = 100.0        # Size of the square domain
const molecule_length = 20.0     # Length of each molecule
```

---

## 📄 License
This project is licensed under the MIT License.

---

## 👤 Author
Created by [Your Name]. Contributions are welcome!
