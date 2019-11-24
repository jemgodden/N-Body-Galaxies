# MPhys-Project-N-Body-Code
The N-Body code written to simulate interacting galaxies, for MPhys project.

Third generation of an N-Body code for modelling galaxy interactions. Each galaxy has a disk of test particles that do not generate a force on other bodies. 

Force on each body is found and used to give that body a new velocity and position, at each time step. Now uses leapfrog method for calculating velocities and positions during a time step.

This code runs a simulation and plots it directly afterwards. Outputs points at image time step and paths to text file to be plotted. Centres images.

Interaction options in IntData.
