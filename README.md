# MPhys-Project-N-Body-Code
The N-Body code written to simulate interacting galaxies, for MPhys project.

Sixth generation of an N-Body code for modelling galaxy interactions. Each galaxy has a disk of test particles that do not generate a force on other bodies. The disk plane can be rotated in the 3d area using normalised spin values in IntData.py.

Force on each body is found and used to give that body a new velocity and position, at each time step. Uses velocity verlet leapfrog method for calculating velocities and positions during a time step.

This code runs a simulation and outputs particle positions at image time step and paths to text file to be plotted.

Interaction options in IntData.

Images plotted by separate function.

Dark matter halo potential for each galaxy implemented, acting on all particles other than galaxy centre.

Added disk plane rotations.

Added stable randomly distributed galaxy disks.

Added Dynamical Friction.

Stripped NBody Code and automation for finding initial conditions.

Density plotters for galaxy and disk analysis.
