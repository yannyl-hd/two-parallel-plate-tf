# Transient Flow Between Two Parallel Plates

This project investigates the transient flow of a fluid initially at rest between two stationary, parallel plates separated by a distance \(2h\). The flow is driven by a constant pressure gradient.

The simplified dimensionless governing equation is

\[
\frac{\partial U}{\partial \tau}
=
1+
\frac{\partial^2 U}{\partial Y^2},
\]

where

\[
Y=\frac{y}{h},
\qquad
\tau=\frac{\nu t}{h^2},
\qquad
U=\frac{u}{U_{\infty}}.
\]

The boundary and initial conditions are

\[
U(-1,\tau)=U(1,\tau)=0,
\qquad
U(Y,0)=0.
\]

## Numerical Methods

Three numerical methods were implemented:

- Explicit Euler
- Implicit Euler
- Fourth-order Runge--Kutta

The numerical solutions were compared with the analytical solution for the transient velocity field.

## Results

The project investigates the effects of:

- Spatial grid spacing, \(\Delta Y\)
- Time-step size, \(\Delta\tau\)
- Diffusion number,

\[
r=\frac{\Delta\tau}{(\Delta Y)^2}.
\]

The following quantities were evaluated:

- Velocity profile
- Dimensionless flow rate
- Dimensionless wall shear
- Maximum absolute error

The explicit Euler method is stable for

\[
r\leq\frac{1}{2}.
\]

Grid refinement generally reduced the numerical error. The RK4 method produced the smallest velocity-profile errors on the finer grids, while the implicit Euler method provided greater stability than the explicit Euler method. Wall-shear errors remained larger because wall shear was calculated from the numerical velocity gradient at the wall.
