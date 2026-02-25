# Nonlinear Control of an Inverted Pendulum System

# Authors
Ali Eidizadeh — ali80ei@gmail.com

# Projects Overview

This project presents a comparative study of four control strategies for stabilizing an inverted pendulum system:

Linear Quadratic Regulator (LQR)

Sliding Mode Control (SMC)

Terminal Sliding Mode Control (TSMC)

Integral Sliding Mode Control (ISMC)

The main objective is to:

Maintain the pendulum in the upright position

Control the cart position

Ensure robustness under varying reference signals and external disturbances

The inverted pendulum is a classic benchmark for nonlinear and unstable systems, making it an ideal platform for evaluating control performance.

# Objectives

Stabilize the inherently unstable inverted pendulum

Achieve accurate reference tracking

Improve transient response performance

Enhance robustness against disturbances

Compare linear and nonlinear control approaches

# Controllers Implemented
1️⃣ LQR (Linear Quadratic Regulator)

A linear optimal control approach designed based on system linearization.

Performance:

Acceptable performance for simple and constant inputs

Works well near the linear operating point

Limited robustness under time-varying references

Reduced performance under strong nonlinear dynamics

Conclusion:
Suitable for nearly linear systems and simple operating conditions, but limited for highly nonlinear and uncertain environments.

2️⃣ SMC (Sliding Mode Control)

A nonlinear robust control technique based on switching control laws.

Performance:

Strong robustness against disturbances

Good stabilization capability

Exhibits chattering phenomenon in the control input

Conclusion:
Highly robust but suffers from input chattering.

3️⃣ TSMC (Terminal Sliding Mode Control)

An enhanced sliding mode method designed for finite-time convergence.

Performance:

Faster convergence compared to classical SMC

Reduced oscillations in control input

Improved transient response

Conclusion:
Provides a good balance between robustness and reduced chattering.

4️⃣ ISMC (Integral Sliding Mode Control)

An advanced sliding mode technique with improved robustness and steady-state performance.

Performance:

Best overall performance among all methods

Shortest settling time

Fast rise time

Near-zero steady-state error

Accurate reference tracking

Strong robustness against disturbances

Significantly reduced chattering

Conclusion:
ISMC demonstrated superior performance in all scenarios and is highly suitable for nonlinear systems with uncertainties.
