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

PriGalCen = [[], [], []]  # Editable list of x, y and z position coordinates for primary galaxy.
SecGalCen = [[], [], []]  # Editable list of x, y and z position coordinates for secondary galaxy.


def path_read():  # Reads information about the path of the galaxies in the interaction.
    print("\nReading files...\n\n")
    if rewind:
        file1 = open("Backwards/RWPriGalPath.txt", "r")
    else:
        file1 = open("Forwards/PriGalPath.txt",
                     "r")  # Opens primary galaxy file containing positions at all time steps.
    for line in file1:  # Goes through file line by line.
        data1 = line.strip().split()
        PriGal[0].append(float(data1[0]))
        PriGal[1].append(float(data1[1]))  # Appends each x, y and z value to list to be plotted later.
        PriGal[2].append(float(data1[2]))
    file1.close()

    if rewind:
        file2 = open("Backwards/RWSecGalPath.txt", "r")
    else:
        file2 = open("Forwards/SecGalPath.txt",
                     "r")  # Opens secondary galaxy file containing positions at all time steps.
    for line in file2:  # Goes through file line by line.
        data2 = line.strip().split()
        SecGal[0].append(float(data2[0]))
        SecGal[1].append(float(data2[1]))  # Appends each value to list to be plotted later.
        SecGal[2].append(float(data2[2]))
    file2.close()

    if centre_mid:
        path_centring_mid()


def path_centring_mid():  # Centres the path in the images of the interaction to a central point of the two galaxies.
    dx = []
    dy = []  # List of the change in position that each particle will undergo to be centred.
    dz = []
    for k in range(len(PriGal[0])):
        gx = (PriGal[0][k] + SecGal[0][k]) / 2
        gy = (PriGal[1][k] + SecGal[1][k]) / 2  # Calculates point equidistant from two galaxies in
        gz = (PriGal[2][k] + SecGal[2][k]) / 2  # in order to centre image in between them.
        dx.append(gx)
        dy.append(gy)  # Making lists of changes in position.
        dz.append(gz)
    for j in range(len(PriGal[0])):
        PriGal[0][j] += - dx[j]
        SecGal[0][j] += - dx[j]
        PriGal[1][j] += - dy[j]
        SecGal[1][j] += - dy[j]  # Altering positions of each particle to centre the images.
        PriGal[2][j] += - dz[j]
        SecGal[2][j] += - dz[j]


def point_centring_mid():  # Function to centre particles in the images about the mid-point of the two galaxies.
    dx = (x[0] + x[1]) / 2
    dy = (y[0] + y[1]) / 2  # Finding the mid-point of the two galaxies.
    dz = (z[0] + z[1]) / 2
    for i in range(len(x)):
        x[i] += - dx
        y[i] += - dy  # Altering the position of each particle in the simulation.
        z[i] += - dz


def centring_pri():  # A function to centre the images of the interaction on the primary galaxy.
    dx = x[0]
    dy = y[0]
    dz = z[0]
    for k in range(len(PriGal[0])):
        PriGalCen[0].append(PriGal[0][k])
        PriGalCen[1].append(PriGal[1][k])
        PriGalCen[2].append(PriGal[2][k])  # Appends galaxy paths to editable list.
        SecGalCen[0].append(SecGal[0][k])
        SecGalCen[1].append(SecGal[1][k])
        SecGalCen[2].append(SecGal[2][k])
    for i in range(len(x)):
        x[i] += - dx
        y[i] += - dy  # Altering the position of each particle in the simulation.
        z[i] += - dz
    for j in range(len(PriGal[0])):
        PriGalCen[0][j] += - dx
        PriGalCen[1][j] += - dy
        PriGalCen[2][j] += - dz  # Altering positions of each galaxy path to centre the images.
        SecGalCen[0][j] += - dx
        SecGalCen[1][j] += - dy
        SecGalCen[2][j] += - dz


def point_read(title):  # Reads information about all particles in the simulation.
    if rewind:
        file = open("Backwards/rimage_%.2f.txt" % title)
    else:
        file = open("Forwards/image_%.2f.txt" % title)  # Reads file containing positions of each particle at a particular time.

    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        x.append(float(data[2]))
        y.append(float(data[3]))  # Appends each x, y and z value to list to be plotted later.
        z.append(float(data[4]))
        colour.append(data[8])  # Appends colour of each particle to be plotted later.

    if centre_mid:
        point_centring_mid()  # For centring on the mid-point of the two galaxies.
    if centre_pri:
        centring_pri()  # For centring on the primary galaxy.


def info():  # Prints information on interaction using files made during simulation.
    total_rp1 = 0
    total_rp2 = 0
    total_part = 0
    pericentre = 0
    time = 0
    if rewind:
        file = open("Backwards/RewindPericentreInfo.txt")
    else:
        file = open("Forwards/PericentreInfo.txt")  # Opens file containing information on the pericentre of the interaction.
    for line in file:  # Goes through file line by line.
        data = line.strip().split()
        total_rp1 = int(data[0])
        total_rp2 = int(data[1])
        total_part = int(data[2])
        pericentre = float(data[3])  # Pericentre of the interaction.
        time = float(data[4])  # Time at which pericentre occurs during interaction.

    print("There are", total_rp1, "particles in the primary galaxy's disk.\n")  # Prints total number of particles in a
                                                                                # galaxy's disk.
    print("There are", total_rp2, "particles in the primary galaxy's disk.\n")  # Prints total number of particles in a
                                                                                # galaxy's disk.
    print("There are a total of", total_part, "particles in the simulation.\n")  # Prints total number of particles
                                                                                    # in a simulation.
    if secondary_gal:
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

        if centre_pri:
            ax.plot3D(PriGalCen[0], PriGalCen[1], PriGalCen[2], path1)  # Prints edited paths for when centring images on the primary galaxy.
            if secondary_gal:
                ax.plot3D(SecGalCen[0], SecGalCen[1], SecGalCen[2], path2)
        else:
            ax.plot3D(PriGal[0], PriGal[1], PriGal[2], path1)  # Plotting path of primary galaxy.
            if secondary_gal:
                ax.plot3D(SecGal[0], SecGal[1], SecGal[2], path2)  # Plotting path of secondary galaxy.

        for j in range(len(x)):  # Plotting all objects on figure.
            ax.plot3D([x[j]], [y[j]], [z[j]], colour[j])

        if centre_pri:
            PriGalCen[0].clear()
            PriGalCen[1].clear()
            PriGalCen[2].clear()  # Clears editable list of galaxy paths, to be used for next image.
            SecGalCen[0].clear()
            SecGalCen[1].clear()
            SecGalCen[2].clear()
        x.clear()
        y.clear()  # Clears all the lists used to plot image, before appending new values for next image.
        z.clear()
        colour.clear()


def main():  # Calling all functions in order.

    if centre_mid and centre_pri:  # Check to make sure both centring's aren't being used at the same time.
        print("\nError. Please select only one centring.")
        return

    if centre_mid and not secondary_gal:
        print("\nError. Cannot centre on the middle of the interactions if there is only one galaxy.")
        return

    path_read()

    info()

    plot()

    plt.show()


if __name__ == '__main__':
    main()
