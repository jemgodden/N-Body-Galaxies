import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

x = []
y = []  # x, y and z coordinates of particle being plotted.
z = []
colour = []  # Colour of particle in plot.

PriGal = [[], [], []]  # List of x, y and z coordinates for position of primary galaxy during interaction.
SecGal = [[], [], []]  # List of x, y and z coordinates for position of secondary galaxy during interaction.


def path_read():  # Reads information about the path of the galaxies in the interaction.
    print("\nReading files...\n\n")
    file1 = open("PriGalPath.txt")  # Opens primary galaxy file containing positions at all time steps.
    for line in file1:  # Goes through file line by line.
        data1 = line.strip().split()
        PriGal[0].append(float(data1[0]))
        PriGal[1].append(float(data1[1]))  # Appends each x, y and z value to list to be plotted later.
        PriGal[2].append(float(data1[2]))

    file2 = open("SecGalPath.txt")  # Opens secondary galaxy file containing positions at all time steps.
    for line in file2:  # Goes through file line by line.
        data2 = line.strip().split()
        SecGal[0].append(float(data2[0]))
        SecGal[1].append(float(data2[1]))  # Appends each value to list to be plotted later.
        SecGal[2].append(float(data2[2]))


def point_read(title):  # Reads information about all particles in the simulation.
    file = open("image_%.2f.txt" % title)  # Reads file containing positions of each particle at a particular time.
    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        x.append(float(data[2]))
        y.append(float(data[3]))  # Appends each x, y and z value to list to be plotted later.
        z.append(float(data[4]))
        colour.append(data[8])  # Appends colour of each particle to be plotted later.
    if centre_mid == True:
        point_centring()


def info():  # Prints information on interaction using files made during simulation.
    pericentre = 0
    time = 0
    print("There are", tot_rp, "particles in one galaxy's disk.\n")  # Prints total number of particles in a
                                                                     # galaxy's disk.
    print("There are a total of", tot_part, "particles in the simulation.\n")  # Prints total number of particles
                                                                               # in a simulation.
    file = open("PericentreInfo.txt")  # Opens file containing information on the pericentre of the interaction.
    for line in file:  # Goes through file line by line.
        data = line.strip().split()
        pericentre = float(data[0])  # Pericentre of the interaction.
        time = float(data[1])  # Time at which pericentre occurs during interaction.

    print("Pericentre of interaction =", round(pericentre, 2), "kpc, occurring at t =", round(time, 2), "Gyrs.\n")
                                                            # Prints value of pericentre and time at which it occurs.
    print("The time between each of the", (frames + 1), "images is", image_time_step, "Gyrs.\n\n")
                                                            # Prints number of images and time between them.


def plot():  # Plot images of interaction.
    print("Producing images...\n")
    for i in range(frames, -1, -1):  # Plots data of each time step on a separate figure.
        title = i * image_time_step  # (frames - i)
        point_read(title)
        fig = plt.figure()  # Creates figure.
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.set_title('t = %.2f Gyrs' % title, fontsize=10)  # Title the figure with the time of interaction.
        ax.set_xlabel('X', fontsize=10)
        ax.set_ylabel('Y', fontsize=10)  # Set axis labels.
        ax.set_zlabel('Z', fontsize=10)
        ax.set_xlim(-2e21, 2e21)
        ax.set_ylim(-2e21, 2e21)  # Set axis limits.
        ax.set_zlim(-2e11, 2e11)
        ax.plot3D(PriGal[0], PriGal[1], PriGal[2], path1)  # Plotting path of primary galaxy.
        ax.plot3D(SecGal[0], SecGal[1], SecGal[2], path2)  # Plotting path of secondary galaxy.
        for j in range(len(x)):  # Plotting all objects on figure.
            ax.plot3D([x[j]], [y[j]], [z[j]], colour[j])
        x.clear()
        y.clear()  # Clears all the lists used to plot image, before appending new values for next image.
        z.clear()
        colour.clear()


def main():  # Calling all functions in order.

    path_read()

    info()

    plot()

    plt.show()


if __name__ == '__main__':
    main()
