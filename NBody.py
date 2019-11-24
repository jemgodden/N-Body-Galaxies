import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

objects = []  # List of all objects in galaxy.


class Body:
    def __init__(self, name, m, x, y, z, vx, vy, vz, colour):
        self.name = name  # Name of body.
        self.m = m  # Mass of body.
        self.x = x  # x position of body.
        self.y = y  # y position of body.
        self.z = z  # z position of body.
        self.vx = vx  # x velocity of body.
        self.vy = vy  # y velocity of body.
        self.vz = vz  # z velocity of body.
        self.saved_x = []  # Array of all x position values of body.
        self.saved_y = []  # Array of all y position values of body.
        self.saved_z = []  # Array of all z position values of body.
        self.image_x = []  # Array of all x position values of body, used for images.
        self.image_y = []  # Array of all y position values of body, used for images.
        self.image_z = []  # Array of all z position values of body, used for images.
        self.colour = colour  # Colour of body on images.

    def name(self):
        return self.name

    def mass(self):
        return self.m

    def x(self):
        return self.x

    def y(self):
        return self.y

    def z(self):
        return self.z

    def vx(self):
        return self.vx

    def vy(self):
        return self.vy

    def vz(self):
        return self.vz

    def saved_x(self):
        return self.saved_x

    def saved_y(self):
        return self.saved_y

    def saved_z(self):
        return self.saved_z

    def image_x(self):
        return self.image_x

    def image_y(self):
        return self.image_y

    def image_z(self):
        return self.image_z

    def colour(self):
        return self.colour

    def force_calc(self, other):  # Force calculation between two bodies
        rx = self.x - other.x
        ry = self.y - other.y  # Distance between two bodies in all directions.
        rz = self.z - other.z
        f = -(G * self.m * other.m)/(rx ** 2 + ry ** 2 + rz ** 2)  # Total force calculation.
        theta = math.atan2(ry, rx)  # Azimuthal angle.
        phi = math.acos(rz / math.sqrt(rx ** 2 + ry ** 2 + rz ** 2))  # Polar angle.
        fx = math.cos(theta) * math.sin(phi) * f
        fy = math.sin(theta) * math.sin(phi) * f  # Force calculation in each direction.
        fz = math.cos(phi) * f
        return fx, fy, fz


def create_galaxies():  # Creates two interacting galaxies using Body class.
    objects.append(Body("Primary", mg1, xg1, yg1, zg1, vxg1, vyg1, vzg1, 'bo'))
    # Creates primary galaxy and adds it to list of bodies.
    objects.append(Body("Secondary", mg2, xg2, yg2, zg2, vxg2, vyg2, vzg2, 'ro'))
    # Creates secondary galaxy and adds it to list of bodies.


def create_rings():  # Creates set of rings for each galaxy.
    for k in range(0, no_rings):  # Creating rings for primary galaxy.
        r = (k + 1) * ring_rad  # Radius of each ring.
        n = (k + 1) * no_rp  # Number of particles in each ring.
        v = math.sqrt(G * mg1 / r)  # Velocity of particles in each ring.
        for i in range(0, n):
            angle = 2 * math.pi * i / n  # Assigning particles to a ring formation.
            x = r * math.cos(angle)
            y = r * math.sin(angle)  # x and y coordinates of each particle in ring.
            vx = v * math.sin(angle)
            vy = -v * math.cos(angle)  # x and y velocities of each particle in ring.
            objects.append(Body("Test", 1, xg1 + x, yg1 + y, zg1, vxg1 + vx, vyg1 + vy, vzg1, 'c.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.

    for k in range(0, no_rings):  # Creating rings for secondary galaxy.
        r = (k + 1) * ring_rad  # Radius of each ring.
        n = (k + 1) * no_rp  # Number of particles in each ring.
        v = math.sqrt(G * mg2 / r)  # Velocity of particles in each ring.
        for i in range(0, n):
            angle = 2 * math.pi * i / n  # Assigning particles to a ring formation.
            x = r * math.cos(angle)
            y = r * math.sin(angle)  # x and y coordinates of each particle in ring.
            vx = v * math.sin(angle)
            vy = -v * math.cos(angle)  # x and y velocities of each particle in ring.
            objects.append(Body('Test', 1, xg2 + x, yg2 + y, zg2, vxg2 + vx, vyg2 + vy, vzg2, 'm.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.


def leapfrog_initial(bodies):  # Produces kick start for the leapfrog algorithm.
    print("Calculating...")
    for body in bodies:
        body.saved_x.append(body.x)
        body.saved_y.append(body.y)  # Appends initial positions of each body to a
        body.saved_z.append(body.z)  # list of saved positions for that body.

    force = {}
    for body in bodies:
        total_fx = total_fy = total_fz = 0.0  # Sets force for all bodies to 0.
        for other in bodies:
            if body is other:  # Checking that not calculating force of a body on itself.
                continue
            if other.name != "Test":  # Does not calculate force due to ring/test particles.
                fx, fy, fz = body.force_calc(other)
                total_fx += fx
                total_fy += fy  # Add together forces of al other bodies acting on that body.
                total_fz += fz
        force[body] = (total_fx, total_fy, total_fz)

    for body in bodies:  # Kick start position and velocities for all bodies.
        fx, fy, fz = force[body]

        body.vx += (fx / body.m) * (time_step / 2)
        body.vy += (fy / body.m) * (time_step / 2)  # Calculates initial half step for velocity in each direction,
        body.vz += (fz / body.m) * (time_step / 2)  # using v_0.5 = v_0 + a * (dt/2), where a = F(x_0)/m.

        body.x += body.vx * time_step
        body.y += body.vy * time_step  # Uses half step in velocity to calculate new position,
        body.z += body.vz * time_step  # using x_1 = x_0 + v_0.5 * dt.

        body.saved_x.append(body.x)
        body.saved_y.append(body.y)  # Saving new position to list of previous positions of body.
        body.saved_z.append(body.z)


def leapfrog(bodies):  # Updates the position and velocity of each particle using a leapfrog algorithm.
    step = 1
    percent = 0.0
    print(percent)
    while True:
        step += 1
        for a in range(100):
            if step == (a * no_step) / 100:
                percent += 1.0
                print(percent)  # Percentage calculator for running the simulation.
        if step == no_step:
            return  # Stop simulation when all steps done.
        else:
            force = {}
            for body in bodies:
                total_fx = total_fy = total_fz = 0.0  # Sets force for all bodies to 0.
                for other in bodies:
                    if body is other:  # Checking that not calculating force of a body on itself.
                        continue
                    if other.name != "Test":  # Does not calculate force due to ring/test particles.
                        fx, fy, fz = body.force_calc(other)
                        total_fx += fx
                        total_fy += fy  # Add together forces of al other bodies acting on that body.
                        total_fz += fz
                force[body] = (total_fx, total_fy, total_fz)

            for body in bodies:
                fx, fy, fz = force[body]

                body.vx += (fx / body.m) * time_step
                body.vy += (fy / body.m) * time_step  # Calculates the new velocity in each direction,
                body.vz += (fz / body.m) * time_step  # using v_i+1.5 = v_i+0.5 + a * dt, where a = F(x_i)/m.

                body.x += body.vx * time_step
                body.y += body.vy * time_step  # Uses new velocity to calculate new position,
                body.z += body.vz * time_step  # using x_i+1 = x_i + v_i+0.5 * dt.

                body.saved_x.append(body.x)
                body.saved_y.append(body.y)  # Saving new position to list of previous positions of body.
                body.saved_z.append(body.z)


def centring(bodies):  # Centres the images of the interaction to a certain point.
    dx = []
    dy = []  # List of the change in position that each particle will undergo to be centred.
    dz = []
    for k in range(len(objects[0].saved_x)):
        x = (objects[0].saved_x[k] + objects[1].saved_x[k]) / 2
        y = (objects[0].saved_y[k] + objects[1].saved_y[k]) / 2  # Calculates point equidistant from tow galaxies in
        z = (objects[0].saved_z[k] + objects[1].saved_z[k]) / 2  # in order to centre image in between them.
        dx.append(x)
        dy.append(y)  # Making lists of changes in position.
        dz.append(z)
    for body in bodies:
        for j in range(len(objects[0].saved_x)):
            body.saved_x[j] += - dx[j]
            body.saved_y[j] += - dy[j]  # Altering positions of each particle to centre the images.
            body.saved_z[j] += - dz[j]
    for body in bodies:
        for i in range(len(objects[0].saved_x)):
            if i % (interval - 1) == 0 or i == no_step:
                body.image_x.append(body.saved_x[i])
                body.image_y.append(body.saved_y[i])  # Saving new position to list of positions that will be plotted.
                body.image_z.append(body.saved_z[i])


def info():
    print("\n\nThere are", tot_rp, "particles in one galaxy's disk.\n")  # Prints total number of particles in a
                                                                         # galaxy's disk.
    print("There are a total of", tot_part, "particles in the simulation.\n")  # Prints total number of particles
                                                                               # in a simulation.


def pericentre_calc():  # Calculates the pericentre of the interaction and the time at which it occurs.
    pericentres = []  # Create list of pericentres at each time step.
    for i in range(len(objects[0].saved_x)):
        xpc = (objects[0].saved_x[i] - objects[1].saved_x[i]) ** 2
        ypc = (objects[0].saved_y[i] - objects[1].saved_y[i]) ** 2
        zpc = (objects[0].saved_z[i] - objects[1].saved_z[i]) ** 2
        pericentre = math.sqrt(xpc + ypc + zpc)  # Calculate pericentre.
        pericentres.append(pericentre)
    min_pericentre = min(pericentres)  # Minimum value of pericentre found in list.
    min_pericentre_pc = min_pericentre / kpc
    t = ((pericentres.index(min(pericentres)) + 1) * time_step) / Gyr  # Min pericentre position in list used
                                                                       # to find time.
    print("Pericentre of interaction =", round(min_pericentre_pc, 2), "kpc, occurring at t =", round(t, 2), "Gyrs.\n")
                                                            # Prints value of pericentre and time at which it occurs.
    print("The time between images is", image_time_step, "Gyrs.\n\n")  # Prints time in between images.

    file = open("PericentreInfo.txt", "w+")
    file.write("{0} {1}".format(min_pericentre_pc, t))  # Prints value of pericentre and time it occurs to a file.
    file.close()


def path_print():
    file1 = open("PriGalPath.txt", "w+")  # Prints coordinates of primary galaxy at every time step to a file.
    for i in range(len(objects[0].saved_x)):
        file1.write("{0} {1} {2}\n".format(objects[0].saved_x[i], objects[0].saved_y[i], objects[0].saved_z[i]))
    file1.close()

    file2 = open("SecGalPath.txt", "w+")  # Prints coordinates of secondary galaxy at every time step to a file.
    for i in range(len(objects[1].saved_x)):
        file2.write("{0} {1} {2}\n".format(objects[1].saved_x[i], objects[1].saved_y[i], objects[1].saved_z[i]))
    file2.close()


def position_print(bodies):  # Printing position of all particles at the time an image is seen.
    for i in range(len(objects[0].image_x)):
        txt_title = i * image_time_step
        file = open("image_%.2f.txt" % txt_title, "w+")  # Opens/creates a file with the name of the time of image.
        for body in bodies:
            file.write("{0} {1} {2} {3}\n".format(body.image_x[i], body.image_y[i], body.image_z[i], body.colour))
        file.close()


def plot():  # Plot images of interaction.
    print("Producing Images...")
    for j in range(len(objects)):  # Shows last image first, so plot in reverse to see in order.
        objects[j].saved_x.reverse()
        objects[j].saved_y.reverse()
        objects[j].saved_z.reverse()
        objects[j].image_x.reverse()
        objects[j].image_y.reverse()
        objects[j].image_z.reverse()
    for i in range(len(objects[0].image_x)):  # Plots data of each time step on a separate figure.
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        title = (-i + frames) * image_time_step
        ax.set_title('t = %.2f Gyrs' % title, fontsize=10)  # Title the figure with the time of interaction.
        ax.set_xlabel('X', fontsize=10)
        ax.set_ylabel('Y', fontsize=10)  # Set axis labels.
        ax.set_zlabel('Z', fontsize=10)
        ax.set_xlim(-2e21, 2e21)
        ax.set_ylim(-2e21, 2e21)  # Set axis limits.
        ax.set_zlim(-2e11, 2e11)
        ax.plot3D(objects[0].saved_x, objects[0].saved_y, objects[0].saved_z, path1)  # Plotting primary galaxy's path.
        ax.plot3D(objects[1].saved_x, objects[1].saved_y, objects[1].saved_z, path2)  # Plotting secondary galaxy's
                                                                                      # path.
        for j in range(len(objects)):  # Plotting all objects on figure.
            ax.plot3D([objects[j].image_x[i]], [objects[j].image_y[i]], [objects[j].image_z[i]], objects[j].colour)


def main():  # Calling all functions in order.

    create_galaxies()

    create_rings()

    leapfrog_initial(objects)

    leapfrog(objects)

    centring(objects)

    info()

    pericentre_calc()

    path_print()

    position_print(objects)

    plot()

    plt.show()


if __name__ == '__main__':
    main()
