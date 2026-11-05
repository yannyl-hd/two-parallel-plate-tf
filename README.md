# Transient Flow Between Two Parallel Plates

This project investigates the transient flow of a fluid initially at rest between two stationary, parallel plates separated by a distance \(2h\). The flow is driven by a constant pressure gradient.

## Governing Equation

The dimensionless governing equation is

$$\frac{\partial U}{\partial \tau} = 1+\frac{\partial^2 U}{\partial Y^2} $$

where

$$ Y=\frac{y}{h},
\qquad
\tau=\frac{\nu t}{h^2},
\qquad
U=\frac{u}{U_{\infty}}. $$

The boundary and initial conditions are

$$ U(-1,\tau)=U(1,\tau)=0 $$

and

$$ U(Y,0)=0 $$

## Analytical Solution

The steady-state solution is

$$ U_s(Y)=\frac{1}{2}(1-Y^2). $$

The transient analytical solution is

$$ U(Y,\tau) =
\frac{1}{2}(1-Y^2)
-
\frac{16}{\pi^3}
\sum_{n=0}^{\infty}
\frac{(-1)^n}{(2n+1)^3}
\cos\left(\frac{(2n+1)\pi Y}{2}\right)
\exp\left[
-\frac{(2n+1)^2\pi^2}{4}\tau
\right]. $$

## Numerical Methods

Three numerical methods were implemented:

- Explicit Euler
- Implicit Euler
- Fourth-order Runge--Kutta (RK4)

The diffusion number is defined as

$$ r=\frac{\Delta\tau}{(\Delta Y)^2}. $$

For the explicit Euler method, the stability condition is

$$ r\leq\frac{1}{2}. $$

## Quantities Evaluated

The numerical solutions were compared with the analytical solution using:

- Dimensionless velocity profiles
- Dimensionless flow rate
- Dimensionless wall shear
- Maximum absolute error

The dimensionless flow rate was calculated by integrating the velocity profile across the channel:

$$ Q(\tau)=\int_{-1}^{1}U(Y,\tau)\,dY. $$

The wall shear was calculated from the velocity gradient at the wall:

$$ \tau_w^*
=
\left.\frac{\partial U}{\partial Y}\right|_{Y=-1}. $$

## Grid Studies

The effects of spatial and temporal discretization were investigated by varying:

- Spatial grid spacing, $$\(\Delta Y\) $$
- Time-step size, $$\(\Delta\tau\) $$
- Diffusion number, $$\(r\)$$

The spatial-grid study used

$$ \Delta Y=0.2,\quad 0.1,\quad 0.05 $$

with a constant diffusion number of \(r=0.25\).

A diffusion-number study was also performed using a constant spatial spacing of $$\(\Delta Y=0.1\)$$.

## Results

All three numerical methods produced velocity profiles that closely agreed with the analytical solution. Reducing the spatial grid spacing generally reduced the numerical error in the velocity profile, flow rate, and wall shear.

The RK4 method produced the smallest velocity-profile errors on the finer grids. The implicit Euler method was generally more accurate than the explicit Euler method among the two Euler schemes and provided greater stability.

The wall-shear error remained larger than the velocity-profile and flow-rate errors because wall shear was calculated from a numerical velocity gradient at the wall. Numerical differentiation amplifies discretization errors, especially when only a small number of points is used.

The explicit Euler method remained stable when

$$ r\leq\frac{1}{2}. $$

Changing the diffusion number did not produce a perfectly monotonic error trend because spatial and temporal discretization errors interact. For the RK4 method, the error was nearly constant over several time-step sizes because the fixed spatial grid limited the overall accuracy.

## Repository Contents

- `tf_2_parallelplate.m` — MATLAB code containing the analytical and numerical solutions
- `AE_601_Midterm_Project_Hannah_Douthit.pdf` — Final project report
- `README.md` — Project description and numerical-method summary

## How to Run

1. Download or clone this repository.
2. Open `tf_2_parallelplate.m` in MATLAB.
3. Select the desired values of $\(\Delta Y\)$, $\(r\)$, and $\(\Delta\tau\).$
4. Run the script.
5. The velocity-profile, flow-rate, wall-shear, and error plots will be generated.

## Author
Hannah Douthit
Hannah Douthit  
AE 601
