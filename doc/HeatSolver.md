---
title: "HeatSolver"
header-includes:
  - \usepackage{amsmath, amsthm, amssymb, xfrac}
  - \usepackage{hyperref}
  - \hypersetup{colorlinks=true}
  - \usepackage[nameinlink,noabbrev]{cleveref}
---

HeatSolver
==========

Solver for two dimensional steady-state heat equation,

\begin{equation}
\label{eq:steady-heat}
 \alpha \nabla^2T = 0 \\
\end{equation}
\begin{equation*}
\alpha = \frac{k}{C_p \rho}
\end{equation*}

where

  - $\alpha$ is the [thermal diffusivity](https://en.wikipedia.org/wiki/Thermal_diffusivity) ($\sfrac{\mathrm m^2}{\mathrm s}$)
  - $k$ is the thermal conductivity ($\sfrac{\mathrm W}{\mathrm m \cdot \mathrm K}$)
  - $C_p$ is the specific heat capacity ($\sfrac{\mathrm J}{\mathrm{kg} \cdot \mathrm K}$)
  - $\rho$ is the density ($\sfrac{\mathrm{kg}}{\mathrm m^3}$)

Solver
======

The implemented solvers are based on the following algebraic equation:

From \cref{eq:steady-heat}:
\begin{align*}
0 &= \frac{\alpha}{\Delta x^2}(T_{i+1,j} + T_{i-1,j} - 2T_{i,j}) + \frac{\alpha}{\Delta y^2}(T_{i+1,j+1} + T_{i,j-1} - 2T_{i,j}) \\
\frac{2T_{i,j}(\Delta x^2 + \Delta y^2)}{\Delta x^2 \Delta y^2} &= \frac{\alpha}{\Delta x^2}(T_{i+1,j} + T_{i-1,j}) + \frac{\alpha}{\Delta y^2}(T_{i+1,j+1} + T_{i,j-1}) \\
T_{i,j} &= \frac{\alpha \Delta y^2}{2(\Delta x^2 + \Delta y^2)}(T_{i+1,j} + T_{i-1,j}) + \frac{\alpha \Delta x^2}{2(\Delta x^2 + \Delta y^2)}(T_{i+1,j+1} + T_{i,j-1})
\end{align*}

In the particular case where $\Delta x = \Delta y$, the equation simplifies and we obtain:

\begin{equation}
\frac{\alpha}{4}(T_{i+1,j} + T_{i-1,j} + T_{i+1,j+1} + T_{i,j-1}) = T_{i,j}
\end{equation}
