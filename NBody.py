import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

objects = []  # List of all objects in galaxy.

start_time = time.time()  # Sets start time in order to find runtime of program.


class Body:
    def __init__(self, name, m, position, velocity, colour):
        self.name = name  # Name of body.
        self.m = m  # Mass of body.
        self.xyz = position  # Array of x, y and z position of body.
        self.vxyz = velocity  # Array of x, y and z velocities of body.
        self.axyz = [0, 0, 0]  # Array of x, y and z acceleration of body.
        self.saved_xyz = [[], [], []]  # Array of all x, y and z position values of body.
        self.colour = colour  # Colour of body on images.

    def name(self):
        return self.name

    def mass(self):
        return self.m

    def xyz(self):
        return self.xyz

    def vxyz(self):
        return self.vxyz

    def axyz(self):
        return self.axyz

    def saved_xyz(self):
        return self.saved_xyz

    def force(self):
        return self.force

    def colour(self):
        return self.colour

    def accel_calc(self, other):  # Force calculation between two bodies
        rx = self.xyz[0] - other.xyz[0]
        ry = self.xyz[1] - other.xyz[1]  # Distance between two bodies in all directions.
        rz = self.xyz[2] - other.xyz[2]
        r = rx ** 2 + ry ** 2 + rz ** 2
        if softening:
            a = - (G * other.m) / (r + (soft_param ** 2))  # Total force calculation.
        else:
            a = - (G * other.m) / r  # Total force calculation.
        theta = math.atan2(ry, rx)  # Azimuthal angle.
        phi = math.acos(rz / math.sqrt(r))  # Polar angle.
        ax = math.cos(theta) * math.sin(phi) * a
        ay = math.sin(theta) * math.sin(phi) * a  # Force calculation in each direction.
        az = math.cos(phi) * a
        return ax, ay, az


def option_checks():
    if not primary_gal and primary_disk:
        print("\nError. There has to a primary galaxy in order to have a primary galaxy disk.")
        exit(1)
    if not primary_gal and primary_dmh:
        print("\nError. There has to a primary galaxy in order to have a primary galaxy dark matter halo.")
        exit(1)
    if not secondary_gal and secondary_disk:
        print("\nError. There has to a secondary galaxy in order to have a secondary galaxy disk.")
        exit(1)
    if not secondary_gal and secondary_dmh:
        print("\nError. There has to a secondary galaxy in order to have a secondary galaxy dark matter halo.")
        exit(1)
    if primary_isolation and secondary_isolation:
        print("\nError. Only one galaxy can be run in isolation at a time.")
        exit(1)
    if (initial_txt and rewind and galaxy_files) or (initial_txt and rewind) or (initial_txt and galaxy_files) \
            or (rewind and galaxy_files):
        print("\nError. Please choose only one way of reading in files.")
        exit(1)


def make_directories():  # Checks directories exist and makes them if not.
    if not os.path.exists('./Forwards'):
        os.makedirs('./Forwards')
    if not os.path.exists('./Backwards'):
        os.makedirs('./Backwards')


def create_galaxies():  # Creates two interacting galaxies using Body class.
    if primary_gal:
        objects.append(Body("Primary", mg1, [xg1, yg1, zg1], [vxg1, vyg1, vzg1], 'bo'))
        # Creates primary galaxy and adds it to list of bodies.
    if secondary_gal:
        objects.append(Body("Secondary", mg2, [xg2, yg2, zg2], [vxg2, vyg2, vzg2], 'ro'))
        # Creates secondary galaxy and adds it to list of bodies.


def create_rings():  # Creates set of rings for each galaxy.
    if primary_disk:
        for k in range(0, no_rings1):  # Creating rings for primary galaxy.
            r = (k + 1) * ring_rad1  # Radius of each ring.
            n = (k + 1) * no_rp1  # Number of particles in each ring.
            if primary_dmh:
                v = math.sqrt((G * mg1 / r) - ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) *
                              (((r / (r + R_s1)) - (math.log(1 + (r / R_s1)))) / r))  # Velocity of particles in ring.
            else:
                v = math.sqrt(G * mg1 / r)

            for i in range(0, n):
                theta = 2 * math.pi * i / n  # Assigning particles to a ring formation.

                xtemp = r * math.cos(theta)
                ytemp = r * math.sin(theta)
                ztemp = 0

                vxtemp = v * math.sin(theta)
                vytemp = -v * math.cos(theta)  # x and y velocities of each particle in ring.
                vztemp = 0

                if norm_spin1[0] == 0 and norm_spin1[1] == 0 and norm_spin1[2] == 1:
                    x = xtemp
                    y = ytemp
                    z = ztemp
                    vx = vxtemp
                    vy = vytemp
                    vz = vztemp

                else:
                    alpha = math.acos(norm_spin1[2])
                    beta = math.asin(norm_spin1[0]/math.sin(alpha))

                    x = xtemp * math.cos(beta) + (ytemp * math.cos(alpha) + ztemp * math.sin(alpha)) * math.sin(beta)
                    y = - xtemp * math.sin(beta) + (ytemp * math.cos(alpha) + ztemp * math.sin(alpha)) * math.cos(beta)
                    z = - ytemp * math.sin(alpha) + ztemp * math.cos(alpha)

                    vx = vxtemp * math.cos(beta) + (vytemp * math.cos(alpha) + vztemp * math.sin(alpha)) * math.sin(beta)
                    vy = - vxtemp * math.sin(beta) + (vytemp * math.cos(alpha) + vztemp * math.sin(alpha)) * math.cos(beta)
                    vz = - vytemp * math.sin(alpha) + vztemp * math.cos(alpha)

                objects.append(Body("pTest", 1, [xg1 + x, yg1 + y, zg1 + z], [vxg1 + vx, vyg1 + vy, vzg1 + vz], 'c.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.

    if secondary_disk:
        for k in range(0, no_rings2):  # Creating rings for secondary galaxy.
            r = (k + 1) * ring_rad2  # Radius of each ring.
            n = (k + 1) * no_rp2  # Number of particles in each ring.
            if secondary_dmh:
                v = math.sqrt((G * mg2 / r) - ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) *
                              (((r / (r + R_s2)) - (math.log(1 + r / R_s2))) / r))  # Velocity of particles in ring.
            else:
                v = math.sqrt(G * mg2 / r)

            for i in range(0, n):
                theta = 2 * math.pi * i / n  # Assigning particles to a ring formation.

                xtemp = r * math.cos(theta)
                ytemp = r * math.sin(theta)
                ztemp = 0

                vxtemp = v * math.sin(theta)
                vytemp = -v * math.cos(theta)  # x and y velocities of each particle in ring.
                vztemp = 0

                if norm_spin2[0] == 0 and norm_spin2[1] == 0 and norm_spin2[2] == 1:
                    x = xtemp
                    y = ytemp
                    z = ztemp
                    vx = vxtemp
                    vy = vytemp
                    vz = vztemp

                else:
                    alpha = math.acos(norm_spin2[2])
                    beta = math.asin(norm_spin2[0]/math.sin(alpha))

                    x = xtemp * math.cos(beta) + (ytemp * math.cos(alpha) + ztemp * math.sin(alpha)) * math.sin(beta)
                    y = - xtemp * math.sin(beta) + (ytemp * math.cos(alpha) + ztemp * math.sin(alpha)) * math.cos(beta)
                    z = - ytemp * math.sin(alpha) + ztemp * math.cos(alpha)

                    vx = vxtemp * math.cos(beta) + (vytemp * math.cos(alpha) + vztemp * math.sin(alpha)) * math.sin(beta)
                    vy = - vxtemp * math.sin(beta) + (vytemp * math.cos(alpha) + vztemp * math.sin(alpha)) * math.cos(beta)
                    vz = - vytemp * math.sin(alpha) + vztemp * math.cos(alpha)

                objects.append(Body("sTest", 1, [xg2 + x, yg2 + y, zg2 + z], [vxg2 + vx, vyg2 + vy, vzg2 + vz], 'm.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.


def rewind_initial():  # Automatically creates the file with initial conditions to run a simulation backwards.
    open("Initial_Conditions.txt", "w+").close()  # Clears contents of file.

    file_name = - time_run / Gyr
    file1 = open("Forwards/image_%.2f.txt" % file_name)  # Reads file at end of last forward simulation.
    lines = file1.readlines()
    file2 = open("Initial_Conditions.txt", "w")
    file2.writelines(lines)  # Writes file full of information for backwards interaction.
    file1.close()
    file2.close()

    read_initial_conditions()


def read_initial_conditions():  # Reads in initial conditions, of all particles, to the simulation form a text file.
    file = open("Initial_Conditions.txt", "r")
    for line in file:  # Reads file line by line.
        data = line.strip().split()
        objects.append(Body(data[0], float(data[1]), [float(data[2]), float(data[3]), float(data[4])],
                            [float(data[5]), float(data[6]), float(data[7])], data[8]))  # Creates bodies using
    file.close()                                                                         # information read from file.


def read_galaxy(file_name):
    file = open(file_name, "r")
    for line in file:
        data = line.strip().split()
        if file_name == "Primary_Galaxy.txt":
            objects.append(Body(data[0], float(data[1]),
                                [xg1 + float(data[2]), yg1 + float(data[3]), zg1 + float(data[4])],
                                [vxg1 + float(data[5]), vyg1 + float(data[6]), vzg1 + float(data[7])], data[8]))
        if file_name == "Secondary_Galaxy.txt":
            objects.append(Body(data[0], float(data[1]),
                                [xg2 + float(data[2]), yg2 + float(data[3]), zg2 + float(data[4])],
                                [vxg2 + float(data[5]), vyg2 + float(data[6]), vzg2 + float(data[7])], data[8]))
    file.close()


def read_galaxy_files():
    read_galaxy("Primary_Galaxy.txt")
    if secondary_gal:
        read_galaxy("Secondary_Galaxy.txt")


def find_galaxy(bodies, galaxy_name):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for body in bodies:
        if body.name == galaxy_name:
            position = objects.index(body)
            break
    return position


def dmh_acceleration(bodies, step_pe, total_pe):
    if primary_dmh:
        pri = find_galaxy(objects, "Primary")
        for body in bodies:
            if body is objects[pri]:
                continue
            else:
                rx = body.xyz[0] - objects[pri].xyz[0]
                ry = body.xyz[1] - objects[pri].xyz[1]  # Distance between two bodies in all directions.
                rz = body.xyz[2] - objects[pri].xyz[2]
                r = math.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
                a = ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) * \
                    (((r / (r + R_s1)) - (math.log(1 + r / R_s1))) / (r ** 2))
                theta = math.atan2(ry, rx)  # Azimuthal angle.
                phi = math.acos(rz / r)  # Polar angle.
                ax = math.cos(theta) * math.sin(phi) * a
                ay = math.sin(theta) * math.sin(phi) * a  # Force calculation in each direction.
                az = math.cos(phi) * a

                body.axyz[0] += ax
                body.axyz[1] += ay  # Add together forces of all other bodies acting on that body.
                body.axyz[2] += az
                if calc_energy:
                    halo1pe = -(G * M_vir1 / r) * (1 / (math.log(1 + c1) - (c1 / (1 + c1)))) * math.log(1 + (r / R_s1))
                    total_pe += halo1pe * body.m / 2

    if secondary_dmh:
        sec = find_galaxy(objects, "Secondary")
        for body in bodies:
            if body is objects[sec]:
                continue
            else:
                rx = body.xyz[0] - objects[sec].xyz[0]
                ry = body.xyz[1] - objects[sec].xyz[1]  # Distance between two bodies in all directions.
                rz = body.xyz[2] - objects[sec].xyz[2]
                r = math.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
                a = ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) * \
                    (((r / (r + R_s2)) - (math.log(1 + r / R_s2))) / (r ** 2))
                theta = math.atan2(ry, rx)  # Azimuthal angle.
                phi = math.acos(rz / r)  # Polar angle.
                ax = math.cos(theta) * math.sin(phi) * a
                ay = math.sin(theta) * math.sin(phi) * a  # Force calculation in each direction.
                az = math.cos(phi) * a

                body.axyz[0] += ax
                body.axyz[1] += ay  # Add together forces of al other bodies acting on that body.
                body.axyz[2] += az
                if calc_energy:
                    halo2pe = -(G * M_vir2 / r) * (1 / (math.log(1 + c2) - (c2 / (1 + c2)))) * math.log(1 + (r / R_s2))
                    total_pe += halo2pe * body.m / 2
    if calc_energy:
        step_pe.append(total_pe)


def find_acceleration(bodies, step_pe, total_pe):
    for body in bodies:
        for i in range(len(body.axyz)):
            body.axyz[i] = 0
        for other in bodies:
            if body is other:  # Checking that not calculating force of a body on itself.
                continue
            if other.name != "Test":  # Does not calculate force due to ring/test particles.
                ax, ay, az = body.accel_calc(other)
                body.axyz[0] += ax
                body.axyz[1] += ay  # Add together forces of al other bodies acting on that body.
                body.axyz[2] += az
                if calc_energy:
                    r = math.sqrt((body.xyz[0]-other.xyz[0]) ** 2 + (body.xyz[1]-other.xyz[1]) ** 2 +
                                  (body.xyz[2]-other.xyz[2]) ** 2)
                    total_pe += - (G * body.m * other.m) / (2*r)

    if primary_dmh or secondary_dmh:
        dmh_acceleration(bodies, step_pe, total_pe)


def leapfrog_initial(bodies, step, step_ke, step_pe):  # Produces kick start for the leapfrog algorithm.
    print("Calculating...")
    position_print(bodies, step)

    total_ke = 0
    total_pe = 0

    for body in bodies:
        body.saved_xyz[0].append(body.xyz[0])
        body.saved_xyz[1].append(body.xyz[1])  # Appends initial positions of each body to a
        body.saved_xyz[2].append(body.xyz[2])  # list of saved positions for that body.
        if calc_energy:
            v = (body.vxyz[0] ** 2) + (body.vxyz[1] ** 2) + (body.vxyz[2] ** 2)
            total_ke += 0.5 * body.m * v
    if calc_energy:
        step_ke.append(total_ke)

    find_acceleration(bodies, step_pe, total_pe)


def leapfrog(bodies):  # Updates the position and velocity of each particle using a leapfrog algorithm.
    step = 0

    step_ke = []
    step_pe = []

    leapfrog_initial(bodies, step, step_ke, step_pe)

    percent = 0.0
    print(percent)
    while True:
        step += 1
        for a in range(100):
            if step == (a * no_step) / 100:
                percent += 1.0
                print(percent)  # Percentage calculator for running the simulation.

        if step == no_step:
            position_print(bodies, step)
            if calc_energy:
                energy_print(step_ke, step_pe)
            return  # Stop simulation when all steps done.
        else:
            total_ke = 0
            total_pe = 0

            for body in bodies:
                for i in range(len(body.axyz)):
                    body.vxyz[i] += body.axyz[i] * (time_step / 2)  # Calculates the new velocity in each direction.
                    body.xyz[i] += body.vxyz[i] * time_step  # Uses new velocity to calculate new position.
                    body.saved_xyz[i].append(body.xyz[i])  # Saving new position to list of previous positions of body.

            find_acceleration(bodies, step_pe, total_pe)

            for body in bodies:
                for i in range(len(body.axyz)):
                    body.vxyz[i] += body.axyz[i] * (time_step / 2)
                if calc_energy:
                    v = (body.vxyz[0] ** 2) + (body.vxyz[1] ** 2) + (body.vxyz[2] ** 2)
                    total_ke += 0.5 * body.m * v
            if calc_energy:
                step_ke.append(total_ke)

            if step % int(interval) == 0:
                position_print(bodies, step)  # Print information on particles to a file at a particular time.


def energy_print(step_ke, step_pe):
    if rewind:
        file1 = open("Backwards/RewindKE.txt", "w+")
    else:
        file1 = open("Forwards/KE.txt", "w+")
    for i in range(len(step_ke)):
        file1.write("{0} \n".format(step_ke[i]))
    file1.close()

    if rewind:
        file2 = open("Backwards/RewindPE.txt", "w+")
    else:
        file2 = open("Forwards/PE.txt", "w+")
    for j in range(len(step_pe)):
        file2.write("{0} \n".format(step_pe[j]))
    file2.close()


def info():
    x = pericentre_calc()
    min_pericentre_pc = x[0]
    t = x[1]

    total_time = (time.time() - start_time) / 60  # Calculates runtime in minutes.

    print("\nThe runtime for this simulation was %.2f minutes.\n\n" % total_time)  # Prints runtime.

    if primary_disk:
        print("There are", tot_rp1, "particles in the primary galaxy's disk.\n")  # Prints total number of particles in a
                                                                              # galaxy's disk.
    if secondary_disk:
        print("There are", tot_rp2, "particles in the secondary galaxy's disk.\n")  # Prints total number of particles
                                                                                    # in a galaxy's disk.
    print("There are a total of", tot_part, "particles in the simulation.\n")  # Prints total number of particles
                                                                               # in a simulation.
    if secondary_gal:
        print("Pericentre of interaction =", round(min_pericentre_pc, 2),
              "kpc, occurring at t =", round(t, 2), "Gyrs.\n")
                # Prints value of pericentre and time at which it occurs.
    print("The time between each of the", (frames + 1), "images is", image_time_step, "Gyrs.\n\n")
                                                            # Prints number of images and time between them.
    if rewind:
        file = open("Backwards/RewindPericentreInfo.txt", "w+")  # Writes information to file for backwards interaction.
    else:
        file = open("Forwards/PericentreInfo.txt", "w+")  # Writes information to file for forwards interaction.
    file.write("{0} {1} {2} {3} {4}".format(tot_rp1, tot_rp2, tot_part, min_pericentre_pc, t))
                                                                    # Prints information values to a file.
    file.close()


def pericentre_calc():  # Calculates the pericentre of the interaction and the time at which it occurs.
    pericentres = []  # Create list of pericentres at each time step.
    xyz_pc = [None] * 3
    pri = find_galaxy(objects, "Primary")
    sec = find_galaxy(objects, "Secondary")
    for i in range(len(objects[pri].saved_xyz[0])):
        for j in range(len(xyz_pc)):
            xyz_pc[j] = (objects[pri].saved_xyz[j][i] - objects[sec].saved_xyz[j][i]) ** 2
        pericentre = math.sqrt(xyz_pc[0] + xyz_pc[1] + xyz_pc[2])  # Calculates pericentre.
        pericentres.append(pericentre)

    min_pericentre = min(pericentres)  # Minimum value of pericentre found in list.
    min_pericentre_pc = min_pericentre / kpc
    t = ((pericentres.index(min(pericentres)) + 1) * time_step) / Gyr  # Min pericentre position in list used
                                                                       # to find time.
    return min_pericentre_pc, t


def path_print():
    if rewind:
        file1 = open("Backwards/RWPriGalPath.txt", "w+")
    else:
        file1 = open("Forwards/PriGalPath.txt", "w+")  # Prints coordinates of primary galaxy at every time step to a file.

    pri = find_galaxy(objects, "Primary")

    for i in range(len(objects[pri].saved_xyz[0])):
        file1.write("{0} {1} {2}\n".format(objects[pri].saved_xyz[0][i], objects[pri].saved_xyz[1][i],
                                           objects[pri].saved_xyz[2][i]))
    file1.close()

    if secondary_gal:
        if rewind:
            file2 = open("Backwards/RWSecGalPath.txt", "w+")
        else:
            file2 = open("Forwards/SecGalPath.txt", "w+")  # Prints coordinates of secondary galaxy at every time step to a file.

        sec = find_galaxy(objects, "Secondary")

        for i in range(len(objects[sec].saved_xyz[0])):
            file2.write("{0} {1} {2}\n".format(objects[sec].saved_xyz[0][i], objects[sec].saved_xyz[1][i],
                                               objects[sec].saved_xyz[2][i]))
        file2.close()


def position_print(bodies, step):  # Printing information of each particle at the time an image is seen.
    txt_title = step * time_step / Gyr
    if rewind:
        file = open("Backwards/rimage_%.3f.txt" % txt_title, "w+")
    else:
        file = open("Forwards/image_%.3f.txt" % txt_title, "w+")  # Opens/creates a file with the name of the time of image.
    for body in bodies:
        file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(body.name, body.m, body.xyz[0], body.xyz[1],
                                            body.xyz[2], body.vxyz[0], body.vxyz[1], body.vxyz[2], body.colour))
    file.close()                    # Writes all information on a particle to one line.


def isolation_print(bodies, file_name):
    file = open(file_name, "w+")
    for body in bodies:
        file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(body.name, body.m, body.xyz[0], body.xyz[1],
                                                                  body.xyz[2], body.vxyz[0], body.vxyz[1], body.vxyz[2],
                                                                  body.colour))
    file.close()


def isolation_output(bodies):
    if primary_isolation:
        create_galaxies()
        create_rings()
        leapfrog(objects)
        info()
        isolation_print(bodies, "Primary_Galaxy.txt")

    if secondary_isolation:
        create_galaxies()
        create_rings()
        leapfrog(objects)
        info()
        isolation_print(bodies, "Secondary_Galaxy.txt")


def main():  # Calling all functions in order.

    option_checks()
    make_directories()

    if primary_isolation or secondary_isolation:
        isolation_output(objects)
        return

    if initial_txt:
        read_initial_conditions()
    if galaxy_files:
        read_galaxy_files()
    if rewind:
        rewind_initial()
    else:
        create_galaxies()
        create_rings()

    leapfrog(objects)
    info()
    path_print()

    print("The data has been printed to text files. Please use the Plotter function to see the images.\n")


if __name__ == '__main__':
    main()
