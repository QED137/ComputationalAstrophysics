# ComputationalAstrophysics
## Numerical Analysis of Fundamental Equations in Astrophysics: Exploring the Formation of Stars

Welcome to the Numerical Analysis repository, where we delve into the computational exploration of fundamental equations in astrophysics. This collection of Python code represents a journey into the intricate processes underlying the formation of stars.
## About Repository

The codes provided here are a product of past coursework in Computational Astrophysics during my master's program. While the algorithms and simulations may have originated from earlier periods, they remain relevant and serve as a valuable resource for understanding the complexities involved in the birth of celestial bodies.
## Acknowledgements
This repository stands as a testament to the interdisciplinary nature of computational astrophysics, where the realms of astrophysics, physics, applied mathematics, and computer science converge. The codes presented here have been refined, it was part of my master course's coursework Computational Astrophysics at University of Tubingen 

Feel free to explore, experiment, and contribute to the ongoing exploration of fundamental astrophysical phenomena. Your engagement is key to advancing our understanding of the mesmerizing processes that govern the formation of stars.

Happy coding and stargazing!

# Kepler's Equation Solver and Celestial Mechanics

This Python script explores solutions to Kepler's equation and demonstrates celestial mechanics calculations. It includes implementations of fixed-point iteration and Newton-Raphson methods to solve Kepler's equation, as well as exercises related to celestial bodies' orbital parameters.

## Contents

1. [Kepler's Equation Solver](#keplers-equation-solver)
2. [Celestial Mechanics](#celestial-mechanics)
    - [Exercise 1: Mercury and Halley](#exercise-1-mercury-and-halley)
    - [Exercise 2: Earth and Mars](#exercise-2-earth-and-mars)

## Kepler's Equation Solver

The script provides two methods for solving Kepler's equation:

- **Fixed-Point Iteration**
- **Newton-Raphson Method**

These methods are applied to calculate the eccentric anomaly, which is crucial in celestial mechanics.

## Celestial Mechanics

### Exercise 1: Mercury and Halley

The script explores the orbital parameters of Mercury and Halley's Comet. It calculates the eccentric anomaly, mean anomaly, and the number of iterations required for convergence using both fixed-point iteration and Newton-Raphson methods.

### Exercise 2: Earth and Mars

The second exercise involves calculating the mutual distance between Earth and Mars over time. It uses Kepler's equation solutions and celestial mechanics principles.

# Physics of Stars

## Lane-Emden Equation Solver and Stellar System Analysis

## Introduction

This repository contains also a numerical solver for the Lane-Emden equation, along with stability analysis tools. The Lane-Emden equation is a differential equation that describes the structure of self-gravitating, spherically symmetric polytropic fluids, such as stars.
## Lane-Emden Equation

The Lane-Emden equation is given by:

$$
\frac{1}{\xi^2} \frac{d}{d\xi} \left( \xi^2 \frac{d\theta}{d\xi} \right) + \theta^n = 0
$$

where $\theta(\xi)$ is a dimensionless function of the radial coordinate ξξ and nn is the polytropic index.

## Numerical SOlver 

The numerical solver implemented in this repository utilizes [insert numerical method here, e.g., Runge-Kutta method] to approximate the solution of the Lane-Emden equation. The solver is designed to handle various polytropic indices and initial conditions.


## Dahlquist's Stability Test

Dahlquist's test assesses the stability of numerical methods for solving ordinary differential equations (ODEs). Developed by Björck and Dahlquist, it analyzes methods using a simple linear test equation:

$$y' = \lambda y $$

The stability is determined by examining the behavior of the stability function $R(h\lambda)$ for different values of $( \mu )$ in the complex plane. A stable method ensures the stability function's absolute value is bounded by 1 for all $( \mu )$ in the left-half complex plane.

Useful for one-step methods like Euler and Runge-Kutta, Dahlquist's test helps choose appropriate methods for solving specific ODEs.




## How to Use

1. Ensure you have Python installed on your machine.
2. Clone this repository: `git clone https://github.com/QED137/ComputationalAstrophysics`
3. Navigate to the repository: `cd ComputationalAstrophysics`
4. Run the script: `python Computational Astrophyiscs -Copy1.ipynb`

Feel free to modify the parameters, explore different celestial bodies, and adapt the code for your own experiments.

## Dependencies

- NumPy
- Matplotlib
- SciPy (for Fresnel integral)



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
