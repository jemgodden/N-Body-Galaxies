import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import random
import sys

# Galaxy options:
rewind = True
newtonian_gravity = False
primary_gal = True  # Option to include the primary galaxy in the interaction.
primary_disk = False  # Option to include a disk of particles as part of the Primary Galaxy.
primary_dmh_potential = True  # Option to include a dark matter halo for the primary galaxy.
primary_live_dmh = False
primary_dynamical_friction = True  # Option for secondary particles to feel dynamical friction from the primary dmh.
secondary_gal = True  # Option to include the secondary galaxy in the interaction.
secondary_disk = False  # Option to include a disk of particles as part of the Secondary Galaxy.
secondary_dmh_potential = True  # Option to include a dark matter halo for the secondary galaxy.
secondary_live_dmh = False
secondary_dynamical_friction = True  # Option  for primary particles to feel dynamical friction from the secondary dmh.

# Constants:
G = 6.67e-11  # Gravitational constant.
km_s = 1e3  # Kilometres per second.
pc = 3.086e16  # Parsec in metres.
kpc = 1e3 * pc  # Kilo-parsec in metres.
sm = 1.989e30  # Solar Mass.
yr = 60 * 60 * 24 * 365  # Year in seconds.
Gyr = 1e9 * yr  # Giga-year in seconds.
critical_density = 136 * sm / (kpc ** 3)  # Critical density of the universe.
separation = float(sys.argv[4])  # Desired separation of the two galaxies.
separation_ratio = separation / 601.23

# Simulation time conditions:
time_step = 1e5 * yr  # Time between each step.
time_run = 1.5 * Gyr  # Total time simulation is run for.
if rewind:
    time_step = - time_step  # Time between each step, running backwards.
    time_run = - time_run  # Total time simulation is run for, backwards.
no_step = time_run / time_step  # Total number of steps in simulation.

# Identifiers for each galaxy, used in the code to specify each.
primary = 1
secondary = 2

# Primary galaxy starting conditions:
pri_galaxy_name = "NGC5257"  # Name of the primary galaxy.
mg1 = 0.358e11 * sm  # Mass of primary galaxy.
xg1 = -211.92 * kpc * separation_ratio  # x position of primary galaxy.
yg1 = 547.54 * kpc * separation_ratio  # y position of primary galaxy.
zg1 = 129.46 * kpc * separation_ratio  # y position of primary galaxy.
vxg1 = float(sys.argv[1])  # x velocity of primary galaxy.
vyg1 = float(sys.argv[2])  # y velocity of primary galaxy.
vzg1 = float(sys.argv[3])  # z velocity of primary galaxy.
pri_galaxy_marker = "bo"  # Colour and size of marker for primary galaxy being plotted.

# Primary galaxy dark matter halo conditions:
M_vir1 = 1.12e11 * sm  # Virial mass of the primary galaxy's dark matter halo.
R_s1 = 4.39 * kpc  # Scale radius of primary galaxy's dark matter halo.
c1 = 6.1  # Concentration of primary galaxy's dark matter halo.
R_vir1 = 26.8 * kpc  # Virial radius of primary galaxy's dark matter halo.
rho_zero1 = 0.05 * sm / (pc ** 3)
V_max1 = 325 * km_s

# Secondary galaxy starting conditions:
sec_galaxy_name = "NGC5258"  # Name of the secondary galaxy.
mg2 = 0.406e11 * sm  # Mass of secondary galaxy.
xg2 = 0 * kpc  # x position of secondary galaxy.
yg2 = 0 * kpc  # y position of secondary galaxy.
zg2 = 0 * kpc  # z position of secondary galaxy.
vxg2 = 0 * km_s  # x velocity of secondary galaxy.
vyg2 = 0 * km_s  # y velocity of secondary galaxy.
vzg2 = 0 * km_s  # z velocity of secondary galaxy.
sec_galaxy_marker = "ro"  # Colour and size of marker for primary galaxy being plotted.

# Secondary galaxy dark matter halo conditions:
M_vir2 = 2.46e11 * sm  # Virial mass of the secondary galaxy's dark matter halo.
R_s2 = 8.30 * kpc  # Scale radius of secondary galaxy's dark matter halo.
c2 = 4.3  # Concentration of secondary galaxy's dark matter halo.
R_vir2 = 35.4 * kpc  # Virial radius of secondary galaxy's dark matter halo.
rho_zero2 = 0.71 * sm / (pc ** 3)
V_max2 = 320 * km_s

objects = []  # List of all objects in galaxy.

start_time = time.time()  # Sets start time in order to find runtime of program.


class Body:
    def __init__(self, name, m, position, velocity, colour):
        self.name = name  # Name of body.
        self.m = m  # Mass of body.
        self.xyz = position  # Array of x, y and z position of body.
        self.v_xyz = velocity  # Array of x, y and z velocities of body.
        self.a_xyz = [0, 0, 0]  # Array of x, y and z acceleration of body.
        self.saved_xyz = [[], [], []]  # Array of all x, y and z position values of body.
        self.saved_v_xyz = [[], [], []]  # Array of all x, y and z velocity values of body.
        self.colour = colour  # Colour of body on images.

    def name(self):
        return self.name

    def mass(self):
        return self.m

    def xyz(self):
        return self.xyz

    def v_xyz(self):
        return self.v_xyz

    def a_xyz(self):
        return self.a_xyz

    def saved_xyz(self):
        return self.saved_xyz

    def saved_v_xyz(self):
        return self.saved_v_xyz

    def colour(self):
        return self.colour

    def find_separation(self, other):
        rx = self.xyz[0] - other.xyz[0]
        ry = self.xyz[1] - other.xyz[1]
        rz = self.xyz[2] - other.xyz[2]
        r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2))
        return [rx, ry, rz], r

    def find_relative_velocity(self, other):
        vx = self.v_xyz[0] - other.v_xyz[0]
        vy = self.v_xyz[1] - other.v_xyz[1]
        vz = self.v_xyz[2] - other.v_xyz[2]
        v = math.sqrt((vx ** 2) + (vy ** 2) + (vz ** 2))
        return [vx, vy, vz], v

    def calculate_newtonian_acceleration(self, other):
        r_xyz, r = self.find_separation(other)

        if softening:
            a = - (G * other.m) / ((r ** 2) + (soft_param ** 2))
        else:
            a = - (G * other.m) / (r ** 2)

        for i in range(len(r_xyz)):
            self.a_xyz[i] += a * (r_xyz[i] / r)

    def calculate_dmh_acceleration(self, other, galaxy_id):
        a = 0
        r_xyz, r = self.find_separation(other)

        if galaxy_id == primary:
            a = ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) * \
                (((r / (r + R_s1)) - (math.log(1 + r / R_s1))) / (r ** 2))
        elif galaxy_id == secondary:
            a = ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) * \
                (((r / (r + R_s2)) - (math.log(1 + r / R_s2))) / (r ** 2))

        for i in range(len(r_xyz)):
            self.a_xyz[i] += a * (r_xyz[i] / r)

    def calculate_dynamical_friction(self, other, causing_galaxy_id):
        r_xyz, r = self.find_separation(other)
        v_xyz, v = self.find_relative_velocity(other)

        density_distribution = 0
        v_dispersion = 0

        epsilon = 28.5 * kpc
        ln_lambda = math.log(r / (1.4 * epsilon))

        if causing_galaxy_id == primary:
            density_distribution = rho_zero1 / ((r / R_s1) * (1 + (r / R_s1)) ** 2)
            v_dispersion = V_max1 * ((1.4393 * (r / R_s1) ** 0.354) / (1 + 1.1756 * (r / R_s1) ** 0.725))

        elif causing_galaxy_id == secondary:
            density_distribution = rho_zero2 / ((r / R_s2) * (1 + (r / R_s2)) ** 2)
            v_dispersion = V_max2 * ((1.4393 * (r / R_s2) ** 0.354) / (1 + 1.1756 * (r / R_s2) ** 0.725))

        X = v / ((2 ** 0.5) * v_dispersion)

        a = - ((4 * math.pi * (G ** 2) * self.m * ln_lambda * density_distribution) / (v ** 2)) * (
                math.erf(X) - (2 * X / (math.pi ** 0.5)) * math.exp(-(X ** 2)))

        for i in range(len(self.a_xyz)):
            self.a_xyz[i] += a * (v_xyz[i] / v)

    def calculate_kinetic_energy(self):
        v = (self.v_xyz[0] ** 2) + (self.v_xyz[1] ** 2) + (self.v_xyz[2] ** 2)

        return 0.5 * self.m * (v ** 2)

    def calculate_potential_energy(self, other):
        r_xyz, r = self.find_separation(other)

        return - (G * self.m * other.m) / (2*r)

    def calculate_dmh_potential_energy(self, galaxy_id):
        halo_pe = 0

        if galaxy_id == primary:
            halo_pe = -(G * M_vir1 / r) * (1 / (math.log(1 + c1) - (c1 / (1 + c1)))) * math.log(1 + (r / R_s1))
        elif galaxy_id == secondary:
            halo_pe = -(G * M_vir2 / r) * (1 / (math.log(1 + c2) - (c2 / (1 + c2)))) * math.log(1 + (r / R_s2))

        return halo_pe * (self.m / 2)


#######################################################################################################################


def find_galaxy(galaxy_name):
    position = 0
    for i in range(len(objects)):
        if objects[i].name == galaxy_name:
            position = i
            break
    return position


#######################################################################################################################


def create_galaxies():
    if primary_gal:
        objects.append(Body(pri_galaxy_name, mg1, [xg1, yg1, zg1], [vxg1, vyg1, vzg1], pri_galaxy_marker))
    if secondary_gal:
        objects.append(Body(sec_galaxy_name, mg2, [xg2, yg2, zg2], [vxg2, vyg2, vzg2], sec_galaxy_marker))


#######################################################################################################################


def find_dynamical_friction(bodies, galaxy_id):
    for body in bodies:
        if galaxy_id == primary:
            other = objects[find_galaxy(pri_galaxy_name)]
            if body.name == sec_galaxy_name:
                body.calculate_dynamical_friction(other, galaxy_id)
            else:
                continue

        if galaxy_id == secondary:
            other = objects[find_galaxy(sec_galaxy_name)]
            if body.name == pri_galaxy_name:
                body.calculate_dynamical_friction(other, galaxy_id)
            else:
                continue


def find_dmh_acceleration(bodies, galaxy_list_position, galaxy_id):
    other = objects[galaxy_list_position]
    for body in bodies:
        if body is other:
            continue
        else:
            body.calculate_dmh_acceleration(other, galaxy_id)


def find_all_dmh_accelerations(bodies):
    if primary_dmh_potential:
        pri = find_galaxy(pri_galaxy_name)
        find_dmh_acceleration(bodies, pri, primary)

    if secondary_dmh_potential:
        sec = find_galaxy(sec_galaxy_name)
        find_dmh_acceleration(bodies, sec, secondary)


def find_newtonian_gravitation(bodies):
    for body in bodies:
        for other in bodies:
            if body is other:
                continue
            else:
                body.calculate_newtonian_acceleration(other)


def find_all_accelerations(bodies):
    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.a_xyz[i] = 0

    if newtonian_gravity:
        find_newtonian_gravitation(bodies)

    if primary_dmh_potential or secondary_dmh_potential:
        find_all_dmh_accelerations(bodies)

    if primary_dynamical_friction:
        find_dynamical_friction(bodies, primary)
    if secondary_dynamical_friction:
        find_dynamical_friction(bodies, secondary)


#######################################################################################################################


def initial_leapfrog_step(bodies):
    for body in bodies:
        for i in range(len(body.saved_xyz)):
            body.saved_xyz[i].append(body.xyz[i])
            if body.name == pri_galaxy_name or body.name == sec_galaxy_name:
                body.saved_v_xyz[i].append(body.v_xyz[i])

    find_all_accelerations(bodies)


def leapfrog_step(bodies):
    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.v_xyz[i] += body.a_xyz[i] * (time_step / 2)
            body.xyz[i] += body.v_xyz[i] * time_step
            body.saved_xyz[i].append(body.xyz[i])

    find_all_accelerations(bodies)

    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.v_xyz[i] += body.a_xyz[i] * (time_step / 2)
            if body.name == pri_galaxy_name or body.name == sec_galaxy_name:
                body.saved_v_xyz[i].append(body.v_xyz[i])


def leapfrog_loop(bodies):
    step = 0

    initial_leapfrog_step(bodies)

    while True:
        step += 1

        if step == no_step:
            return

        else:
            leapfrog_step(bodies)


#######################################################################################################################


def find_separations(separations):
    xyz_pc = [None] * 3

    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    for i in range(len(objects[pri].saved_xyz[0])):
        for j in range(len(xyz_pc)):
            xyz_pc[j] = (objects[pri].saved_xyz[j][i] - objects[sec].saved_xyz[j][i]) ** 2
        separations.append(math.sqrt(xyz_pc[0] + xyz_pc[1] + xyz_pc[2]))


def calculate_separation_info():
    separations = []

    find_separations(separations)

    pericentre = min(separations)
    apocentre = max(separations)
    end_sep = separations[-1]
    max_end = end_sep/apocentre
    pericentre = pericentre / kpc

    time_of_pericentre = ((separations.index(min(separations)) + 1) * time_step) / Gyr

    return pericentre, time_of_pericentre, max_end


def return_info():
    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    pericentre, time_of_pericentre, max_end = calculate_separation_info()

    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}".format(
          objects[pri].xyz[0], objects[pri].xyz[1], objects[pri].xyz[2], objects[pri].v_xyz[0], objects[pri].v_xyz[1],
          objects[pri].v_xyz[2], objects[sec].xyz[0], objects[sec].xyz[1], objects[sec].xyz[2], objects[sec].v_xyz[0],
          objects[sec].v_xyz[1], objects[sec].v_xyz[2], pericentre, time_of_pericentre, separation, vxg1, vyg1, vzg1,
          max_end))


#######################################################################################################################


def generate_simulation():
    create_galaxies()

    leapfrog_loop(objects)

    return_info()


#######################################################################################################################


def main():
    generate_simulation()


if __name__ == '__main__':
    main()
