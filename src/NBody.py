import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import random
from mpl_toolkits.mplot3d import Axes3D
from Config import *

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

        # epsilon = (0.98 * (3002 ** -0.26)) * kpc
        epsilon = 28.5 * kpc
        ln_lambda = math.log(r / (1.4 * epsilon))

        if causing_galaxy_id == primary:
            density_distribution = (102 * critical_density) / ((r / R_s1) * (1 + (r / R_s1)) ** 2)
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

        return - (G * self.m * other.m) / (2 * r)

    def calculate_dmh_potential_energy(self, galaxy_id):
        halo_pe = 0

        if galaxy_id == primary:
            halo_pe = -(G * M_vir1 / r) * (1 / (math.log(1 + c1) - (c1 / (1 + c1)))) * math.log(1 + (r / R_s1))
        elif galaxy_id == secondary:
            halo_pe = -(G * M_vir2 / r) * (1 / (math.log(1 + c2) - (c2 / (1 + c2)))) * math.log(1 + (r / R_s2))

        return halo_pe * (self.m / 2)


#######################################################################################################################


def option_checks():
    if (primary_isolation or secondary_isolation) and gal_sep_plot:
        print("\nError. Cannot plot galaxy separation when simulating in isolation.")
        exit(1)
    if not primary_gal and primary_disk:
        print("\nError. There has to a primary galaxy in order to have a primary galaxy disk.")
        exit(1)
    if not primary_gal and primary_dmh_potential:
        print("\nError. There has to a primary galaxy in order to have a primary galaxy dark matter halo.")
        exit(1)
    if not secondary_gal and secondary_disk:
        print("\nError. There has to a secondary galaxy in order to have a secondary galaxy disk.")
        exit(1)
    if not secondary_gal and secondary_dmh_potential:
        print("\nError. There has to a secondary galaxy in order to have a secondary galaxy dark matter halo.")
        exit(1)
    if primary_isolation and secondary_isolation:
        print("\nError. Only one galaxy can be run in isolation at a time.")
        exit(1)
    if (initial_txt and rewind and galaxy_files) or (initial_txt and rewind) or (initial_txt and galaxy_files) \
            or (rewind and galaxy_files):
        print("\nError. Please choose only one way of reading in files.")
        exit(1)


def make_directories():
    if not os.path.exists('../Forwards'):
        os.makedirs('../Forwards')
    if not os.path.exists('../Backwards'):
        os.makedirs('../Backwards')


#######################################################################################################################


def read_file(file_name):
    file = open(file_name, "r")
    for line in file:
        data = line.strip().split()
        objects.append(Body(data[0], float(data[1]), [float(data[2]), float(data[3]), float(data[4])],
                            [float(data[5]), float(data[6]), float(data[7])], data[8]))
    file.close()


def read_galaxy_file_line(file_name, line):
    data = line.strip().split()
    if file_name == "Primary_Galaxy.txt":
        objects.append(Body(data[0], float(data[1]),
                            [xg1 + float(data[2]), yg1 + float(data[3]), zg1 + float(data[4])],
                            [vxg1 + float(data[5]), vyg1 + float(data[6]), vzg1 + float(data[7])], data[8]))
    if file_name == "Secondary_Galaxy.txt":
        objects.append(Body(data[0], float(data[1]),
                            [xg2 + float(data[2]), yg2 + float(data[3]), zg2 + float(data[4])],
                            [vxg2 + float(data[5]), vyg2 + float(data[6]), vzg2 + float(data[7])], data[8]))


def read_galaxy_file(file_name):
    file = open(file_name, "r")
    for line in file:
        read_galaxy_file_line(file_name, line)
    file.close()


def read_galaxy_files():
    if primary_gal:
        read_galaxy("Primary_Galaxy.txt")
    if secondary_gal:
        read_galaxy("Secondary_Galaxy.txt")


def find_galaxy(galaxy_name):
    position = 0
    for i in range(len(objects)):
        if objects[i].name == galaxy_name:
            position = i
            break
    return position


#######################################################################################################################


def read_dmh_file(file_name, galaxy_id):
    file = open(file_name, "r")
    for line in file:
        data = line.strip().split()
        if galaxy_id == primary:
            objects.append(Body(pri_dmh_name, float(data[1]),
                                [xg1 + float(data[2]), yg1 + float(data[3]), zg1 + float(data[4])],
                                [vxg1 + float(data[5]), vyg1 + float(data[6]), vzg1 + float(data[7])], ''))
        elif galaxy_id == secondary:
            objects.append(Body(sec_dmh_name, float(data[1]),
                                [xg2 + float(data[2]), yg2 + float(data[3]), zg2 + float(data[4])],
                                [vxg2 + float(data[5]), vyg2 + float(data[6]), vzg2 + float(data[7])], ''))
    file.close()


def read_dark_matter_halos():
    if primary_live_dmh:
        read_dmh_file(pri_dmh_file, primary)
    if secondary_live_dmh:
        read_dmh_file(sec_dmh_file, secondary)


#######################################################################################################################


def create_galaxies():
    if primary_gal:
        objects.append(Body(pri_galaxy_name, mg1, [xg1, yg1, zg1], [vxg1, vyg1, vzg1], pri_galaxy_marker))
    if secondary_gal:
        objects.append(Body(sec_galaxy_name, mg2, [xg2, yg2, zg2], [vxg2, vyg2, vzg2], sec_galaxy_marker))


def rotate_disk_particle(norm_spin, xyz_temp, v_xyz_temp):
    xyz = []
    v_xyz = []

    alpha = math.acos(norm_spin[2])
    beta = math.asin(norm_spin[0] / math.sin(alpha))

    xyz.append(xyz_temp[0] * math.cos(beta) + (xyz_temp[1] * math.cos(alpha) +
                                               xyz_temp[2] * math.sin(alpha)) * math.sin(beta))
    xyz.append(- xyz_temp[0] * math.sin(beta) + (xyz_temp[1] * math.cos(alpha) +
                                                 xyz_temp[2] * math.sin(alpha)) * math.cos(beta))
    xyz.append(- xyz_temp[1] * math.sin(alpha) + xyz_temp[2] * math.cos(alpha))

    v_xyz.append(v_xyz_temp[0] * math.cos(beta) + (v_xyz_temp[1] * math.cos(alpha) +
                                                   v_xyz_temp[2] * math.sin(alpha)) * math.sin(beta))
    v_xyz.append(- v_xyz_temp[0] * math.sin(beta) + (v_xyz_temp[1] * math.cos(alpha) +
                                                     v_xyz_temp[2] * math.sin(alpha)) * math.cos(beta))
    v_xyz.append(- v_xyz_temp[1] * math.sin(alpha) + v_xyz_temp[2] * math.cos(alpha))

    return xyz, v_xyz


def place_random_disk_particle(galaxy_id, norm_spin, theta, r, v):
    xyz = []
    v_xyz = []

    xyz_temp = [(r * math.cos(theta)), (r * math.sin(theta)), 0]

    v_xyz_temp = [(v * math.sin(theta)), (-v * math.cos(theta)), 0]

    if norm_spin[0] == 0 and norm_spin[1] == 0 and norm_spin[2] == 1:
        for i in range(len(xyz_temp)):
            xyz.append(xyz_temp[i])
            v_xyz.append(v_xyz_temp[i])

    else:
        xyz, v_xyz = rotate_disk_particle(norm_spin, xyz_temp, v_xyz_temp)

    if galaxy_id == primary:
        objects.append(Body(pri_disk_name, mdp1, [xg1 + xyz[0], yg1 + xyz[1], zg1 + xyz[2]],
                            [vxg1 + v_xyz[0], vyg1 + v_xyz[1], vzg1 + v_xyz[2]], pri_disk_marker))

    if galaxy_id == secondary:
        objects.append(Body(sec_disk_name, mdp2, [xg2 + xyz[0], yg2 + xyz[1], zg2 + xyz[2]],
                            [vxg2 + v_xyz[0], vyg2 + v_xyz[1], vzg2 + v_xyz[2]], sec_disk_marker))


def return_random_disk_particle_velocity(galaxy_id, count, r):
    v = 0

    if galaxy_id == primary:
        enclosed_mass1 = mg1 + (count * mdp1)

        if primary_dmh_potential:
            v = math.sqrt((G * enclosed_mass1 / r) - ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) *
                          (((r / (r + R_s1)) - (math.log(1 + (r / R_s1)))) / r))
        else:
            v = math.sqrt(G * enclosed_mass1 / r)

    elif galaxy_id == secondary:
        enclosed_mass2 = mg2 + (count * mdp2)

        if secondary_dmh_potential:
            v = math.sqrt((G * enclosed_mass2 / r) - ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) *
                          (((r / (r + R_s2)) - (math.log(1 + (r / R_s2)))) / r))
        else:
            v = math.sqrt(G * enclosed_mass2 / r)

    return v


def find_random_disk_particle_velocity(galaxy_id, saved_particles):
    for i in range(len(saved_particles[1])):
        count = 0
        for j in range(len(saved_particles[1])):
            if i == j:
                continue
            elif saved_particles[1][j] < saved_particles[1][i]:
                count += 1
        saved_particles[2].append(return_random_disk_particle_velocity(galaxy_id, count, saved_particles[1][i]))


def generate_random_disk_particle_radius(galaxy_id, theta, saved_particles):
    r = 0
    prob = 0

    if galaxy_id == primary:
        r = dr1 * random.uniform(0.05, 1.0)
        prob = (1 / (sigma1 * math.sqrt(2 * math.pi))) * math.exp(-0.5 * (((r / dr1) - mu1) / sigma1) ** 2)

    elif galaxy_id == secondary:
        r = dr2 * random.uniform(0.05, 1.0)
        prob = (1 / (sigma2 * math.sqrt(2 * math.pi))) * math.exp(-0.5 * (((r / dr2) - mu2) / sigma2) ** 2)

    rand_no = random.uniform(0, 1.0)

    if rand_no <= prob:
        saved_particles[1].append(r)

    else:
        generate_random_disk_particle_radius(galaxy_id, theta, saved_particles)


def create_random_disk(galaxy_id):
    saved_particles = [[], [], []]

    if galaxy_id == primary:
        for i in range(0, tot_dp1):
            theta = 2 * math.pi * i / tot_dp1
            saved_particles[0].append(theta)

            generate_random_disk_particle_radius(galaxy_id, theta, saved_particles)

        find_random_disk_particle_velocity(galaxy_id, saved_particles)

        for j in range(len(saved_particles[0])):
            place_random_disk_particle(galaxy_id, norm_spin1, saved_particles[0][j], saved_particles[1][j],
                                       saved_particles[2][j])

    elif galaxy_id == secondary:
        for i in range(0, tot_dp2):
            theta = 2 * math.pi * i / tot_dp2
            saved_particles[0].append(theta)

            generate_random_disk_particle_radius(galaxy_id, theta, saved_particles)

        find_random_disk_particle_velocity(galaxy_id, saved_particles)

        for j in range(len(saved_particles[0])):
            place_random_disk_particle(galaxy_id, norm_spin2, saved_particles[0][j], saved_particles[1][j],
                                       saved_particles[2][j])


def find_ring_disk_particle_velocity(galaxy_id, r):
    v = 0

    if galaxy_id == primary:
        if primary_dmh_potential:
            v = math.sqrt((G * mg1 / r) - ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) *
                          (((r / (r + R_s1)) - (math.log(1 + (r / R_s1)))) / r))
        else:
            v = math.sqrt(G * mg1 / r)

    elif galaxy_id == secondary:
        if secondary_dmh_potential:
            v = math.sqrt((G * mg2 / r) - ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) *
                          (((r / (r + R_s2)) - (math.log(1 + (r / R_s2)))) / r))
        else:
            v = math.sqrt(G * mg2 / r)

    return v


def place_ring_disk_particle(galaxy_id, norm_spin, theta, r):
    xyz = []
    v_xyz = []

    v = find_ring_disk_particle_velocity(galaxy_id, r)

    xyz_temp = [(r * math.cos(theta)), (r * math.sin(theta)), 0]

    v_xyz_temp = [(v * math.sin(theta)), (-v * math.cos(theta)), 0]

    if norm_spin[0] == 0 and norm_spin[1] == 0 and norm_spin[2] == 1:
        for i in range(len(xyz_temp)):
            xyz.append(xyz_temp[i])
            v_xyz.append(v_xyz_temp[i])

    else:
        xyz, v_xyz = rotate_disk_particle(norm_spin, xyz_temp, v_xyz_temp)

    if galaxy_id == primary:
        objects.append(Body(pri_disk_name, mdp1, [xg1 + xyz[0], yg1 + xyz[1], zg1 + xyz[2]],
                            [vxg1 + v_xyz[0], vyg1 + v_xyz[1], vzg1 + v_xyz[2]], pri_disk_marker))

    if galaxy_id == secondary:
        objects.append(Body(sec_disk_name, mdp2, [xg2 + xyz[0], yg2 + xyz[1], zg2 + xyz[2]],
                            [vxg2 + v_xyz[0], vyg2 + v_xyz[1], vzg2 + v_xyz[2]], sec_disk_marker))


def create_ring_disk(galaxy_id):
    if galaxy_id == primary:
        for k in range(0, no_rings1):
            r = (k + 1) * ring_rad1
            n = (k + 1) * no_rp1

            for i in range(0, n):
                theta = 2 * math.pi * i / n
                place_ring_disk_particle(galaxy_id, norm_spin1, theta, r)

    elif galaxy_id == secondary:
        for k in range(0, no_rings2):
            r = (k + 1) * ring_rad2
            n = (k + 1) * no_rp2

            for i in range(0, n):
                theta = 2 * math.pi * i / n
                place_ring_disk_particle(galaxy_id, norm_spin2, theta, r)


def create_galaxy_disks():
    if primary_disk:
        if random_disks:
            create_random_disk(primary)
        else:
            create_ring_disk(primary)

    if secondary_disk:
        if random_disks:
            create_random_disk(secondary)
        else:
            create_ring_disk(secondary)


#######################################################################################################################


def find_dynamical_friction(bodies, galaxy_id):
    for body in bodies:
        if galaxy_id == primary:
            other = objects[find_galaxy(pri_galaxy_name)]
            if body.name == sec_galaxy_name or body.name == sec_disk_name:
                body.calculate_dynamical_friction(other, galaxy_id)
            else:
                continue

        if galaxy_id == secondary:
            other = objects[find_galaxy(sec_galaxy_name)]
            if body.name == pri_galaxy_name or body.name == pri_disk_name:
                body.calculate_dynamical_friction(other, galaxy_id)
            else:
                continue


def find_dmh_acceleration(bodies, galaxy_list_position, galaxy_id, total_energies):
    other = objects[galaxy_list_position]
    for body in bodies:
        if body is other:
            continue
        else:
            body.calculate_dmh_acceleration(other, galaxy_id)
            if calc_energy:
                total_energies[1] += body.calculate_dmh_potential_energy(other)


def find_all_dmh_accelerations(bodies, total_energies):
    if primary_dmh_potential:
        pri = find_galaxy(pri_galaxy_name)
        find_dmh_acceleration(bodies, pri, primary, total_energies)

    if secondary_dmh_potential:
        sec = find_galaxy(sec_galaxy_name)
        find_dmh_acceleration(bodies, sec, secondary, total_energies)


def find_newtonian_gravitation(bodies, total_energies):
    for body in bodies:
        for other in bodies:
            if body is other:
                continue
            else:
                if (other.name != pri_disk_name) and (other.name != sec_disk_name):
                    body.calculate_newtonian_acceleration(other)
                if calc_energy:
                    total_energies[1] += body.calculate_potential_energy(other)
                else:
                    continue


def find_all_accelerations(bodies, total_energies):
    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.a_xyz[i] = 0

    if newtonian_gravity:
        find_newtonian_gravitation(bodies, total_energies)

    if primary_dmh_potential or secondary_dmh_potential:
        find_all_dmh_accelerations(bodies, total_energies)

    if primary_dynamical_friction:
        find_dynamical_friction(bodies, primary)
    if secondary_dynamical_friction:
        find_dynamical_friction(bodies, secondary)


#######################################################################################################################


def file_print_energies(step_ke, step_pe):
    if rewind:
        file = open("Backwards/RewindEnergies.txt", "w+")
    else:
        file = open("Forwards/Energies.txt", "w+")
    for i in range(len(step_ke)):
        file.write("{0} {1}\n".format(step_ke[i], step_pe[i]))
    file.close()


def append_energies(step_ke, step_pe, total_energies):
    step_ke.append(total_energies[0])
    step_pe.append(total_energies[1])


#######################################################################################################################


def file_print_all_particles(bodies, file_name):
    file = open(file_name, "w+")
    for body in bodies:
        file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(body.name, body.m, body.xyz[0], body.xyz[1],
                                                                  body.xyz[2], body.v_xyz[0], body.v_xyz[1],
                                                                  body.v_xyz[2], body.colour))
    file.close()


def time_file_print_particles(bodies, step):
    txt_title = step * time_step / Gyr

    if rewind:
        file_print_all_particles(bodies, "Backwards/rimage_%.5f.txt" % txt_title)
    else:
        file_print_all_particles(bodies, "Forwards/image_%.5f.txt" % txt_title)


#######################################################################################################################


def initial_leapfrog_step(bodies, step, step_ke, step_pe):
    time_file_print_particles(bodies, step)

    total_energies = [0, 0]

    for body in bodies:
        for i in range(len(body.saved_xyz)):
            body.saved_xyz[i].append(body.xyz[i])
            if body.name == pri_galaxy_name or body.name == sec_galaxy_name:
                body.saved_v_xyz[i].append(body.v_xyz[i])

        if calc_energy:
            total_energies[0] += body.calculate_kinetic_energy()

    find_all_accelerations(bodies, total_energies)

    if calc_energy:
        append_energies(step_ke, step_pe, total_energies)


def leapfrog_step(bodies, step, step_ke, step_pe):
    total_energies = [0, 0]

    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.v_xyz[i] += body.a_xyz[i] * (time_step / 2)
            body.xyz[i] += body.v_xyz[i] * time_step
            body.saved_xyz[i].append(body.xyz[i])

    find_all_accelerations(bodies, total_energies)

    for body in bodies:
        for i in range(len(body.a_xyz)):
            body.v_xyz[i] += body.a_xyz[i] * (time_step / 2)
            if body.name == pri_galaxy_name or body.name == sec_galaxy_name:
                body.saved_v_xyz[i].append(body.v_xyz[i])

        if calc_energy:
            total_energies[0] += body.calculate_kinetic_energy()

    if calc_energy:
        append_energies(step_ke, step_pe, total_energies)

    if step % int(interval) == 0:
        time_file_print_particles(bodies, step)


def leapfrog_loop(bodies):
    step = 0
    step_ke = []
    step_pe = []
    percent_time_start = 0

    print("Calculating...")
    initial_leapfrog_step(bodies, step, step_ke, step_pe)

    percent = 0.0
    print(percent)
    while True:
        step += 1
        for a in range(100):
            if step == (a * no_step) / 100:
                percent += 1.0
                print(percent)
                if percent == 1.0:
                    percent_time_start = time.time()
                if percent == 2.0:
                    percent_time_end = time.time()
                    average_percent_time = percent_time_end - percent_time_start
                    print("One percent took ", average_percent_time, " seconds. It is expected for this simulation to "
                          "take a further", (average_percent_time * 98), " seconds. This is equivalent to ",
                          (average_percent_time * 98) / 60, " minutes.")

        if step == no_step:
            time_file_print_particles(bodies, step)
            if calc_energy:
                file_print_energies(step_ke, step_pe)
            return

        else:
            leapfrog_step(bodies, step, step_ke, step_pe)


#######################################################################################################################


def file_print_path(file_name, galaxy_list_position):
    file = open(file_name, "w+")

    for i in range(len(objects[galaxy_list_position].saved_xyz[0])):
        file.write("{0} {1} {2}\n".format(objects[galaxy_list_position].saved_xyz[0][i],
                                          objects[galaxy_list_position].saved_xyz[1][i],
                                          objects[galaxy_list_position].saved_xyz[2][i]))
    file.close()


def file_print_galaxy_paths():
    if primary_gal:
        pri = find_galaxy(pri_galaxy_name)
        if rewind:
            file_print_path("Backwards/RWPriGalPath.txt", pri)
        else:
            file_print_path("../Forwards/PriGalPath.txt", pri)

    if secondary_gal:
        sec = find_galaxy(sec_galaxy_name)
        if rewind:
            file_print_path("Backwards/RWSecGalPath.txt", sec)
        else:
            file_print_path("../Forwards/SecGalPath.txt", sec)


#######################################################################################################################


def plot_galaxy_separation(separations, rel_velocities):
    time_steps = []

    for k in range(len(separations)):
        separations[k] = separations[k] / kpc
        rel_velocities[k] = rel_velocities[k] / km_s
        time_steps.append(k * time_step / Gyr)

    fig, ax = plt.subplots()
    ax.plot(time_steps, separations, 'r-', label='Relative Distance', linewidth=8)
    # ax.set_ylim([0, 300])
    ax.set_xlabel("Time $(Gyrs)$", fontsize=28, weight='bold')
    ax.set_ylabel("Relative Distance $(kpc)$", fontsize=30, color='red', weight='bold')
    ax2 = ax.twinx()
    ax2.plot(time_steps, rel_velocities, 'b-', label='Relative Velocity', linewidth=8)
    # ax2.set_ylim([0, 300])
    ax2.set_ylabel("Relative Velocity $(kms^{-1})$", fontsize=30, color='blue', weight='bold')

    ax.tick_params(labelsize=28, axis='y', colors='red')
    ax.tick_params(labelsize=28, axis='x')
    ax2.tick_params(labelsize=28, colors='blue')

    plt.show()


def find_separations_and_relative_velocity(separations, rel_velocities):
    xyz_pc = [None] * 3
    v_xyz_pc = [None] * 3

    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    for i in range(len(objects[pri].saved_xyz[0])):
        for j in range(len(xyz_pc)):
            xyz_pc[j] = (objects[pri].saved_xyz[j][i] - objects[sec].saved_xyz[j][i]) ** 2
            v_xyz_pc[j] = (objects[pri].saved_v_xyz[j][i] - objects[sec].saved_v_xyz[j][i]) ** 2
        separations.append(math.sqrt(xyz_pc[0] + xyz_pc[1] + xyz_pc[2]))
        rel_velocities.append(math.sqrt(v_xyz_pc[0] + v_xyz_pc[1] + v_xyz_pc[2]))


def calculate_separation_info():
    separations = []
    post_pc_separations = []
    rel_velocities = []

    find_separations_and_relative_velocity(separations, rel_velocities)

    pericentre = min(separations)
    pericentre_position = separations.index(min(separations))
    pericentre = pericentre / kpc

    for k in range(len(separations)):
        if k >= pericentre_position:
            post_pc_separations.append(separations[k])

    apocentre = max(post_pc_separations) / kpc
    print("\n\nThe apocentre was:", apocentre, "kpc.\n")

    time_of_pericentre = ((separations.index(min(separations)) + 1) * time_step) / Gyr
    if gal_sep_plot:
        plot_galaxy_separation(separations, rel_velocities)

    return pericentre, time_of_pericentre


def print_interaction_info():
    pericentre = 0
    time_of_pericentre = 0
    if not primary_isolation and not secondary_isolation:
        pericentre, time_of_pericentre = calculate_separation_info()
    total_time = (time.time() - start_time) / 60

    print("\nRuntime: %.2f minutes.\n" % total_time)
    if primary_disk and not secondary_isolation:
        print("Primary galaxy disk particles: ", tot_dp1, ".\n")
    if secondary_disk and not primary_isolation:
        print("Secondary galaxy disk particles: ", tot_dp2, ".\n")
    print("Total simulation particles: ", tot_part, ".\n")
    if primary_gal and secondary_gal:
        print("Pericentre: ", round(pericentre, 2), "kpc.\n")
        print("Time of pericentre: ", round(time_of_pericentre, 2), "Gyrs.\n")
    print("Number of images produced: ", (frames + 1), ".\n")
    print("Time between images: ", image_time_step, "Gyrs.\n\n")

    if rewind:
        file = open("Backwards/RewindPericentreInfo.txt", "w+")
    else:
        file = open("../Forwards/PericentreInfo.txt", "w+")
    file.write("{0} {1} {2} {3} {4}".format(tot_dp1, tot_dp2, tot_part, pericentre, time_of_pericentre))
    file.close()


#######################################################################################################################


def isolation_simulation():
    create_galaxies()
    create_galaxy_disks()

    leapfrog_loop(objects)
    print_interaction_info()

    if primary_isolation:
        file_print_all_particles(objects, "../data/Primary_Galaxy.txt")
        print("The primary galaxy was simulated in isolation.\n")
    elif secondary_isolation:
        file_print_all_particles(objects, "../data/Secondary_Galaxy.txt")
        print("The secondary galaxy was simulated in isolation.\n")


def file_simulation():
    read_file("Initial_Conditions.txt")

    leapfrog_loop(objects)

    file_print_galaxy_paths()
    print_interaction_info()


def galaxy_files_simulation():
    read_galaxy_file("../data/Primary_Galaxy.txt")
    read_galaxy_file("../data/Secondary_Galaxy.txt")

    leapfrog_loop(objects)

    file_print_galaxy_paths()
    print_interaction_info()


def generate_simulation():
    create_galaxies()
    create_galaxy_disks()

    leapfrog_loop(objects)

    file_print_galaxy_paths()
    print_interaction_info()


#######################################################################################################################


def main():

    option_checks()
    make_directories()

    if primary_isolation or secondary_isolation:
        isolation_simulation()

    elif initial_txt:
        file_simulation()

    elif galaxy_files:
        galaxy_files_simulation()

    else:
        generate_simulation()

    print("The data has been printed to text files. Please use the Plotter function to see the images.\n")


if __name__ == '__main__':
    main()
