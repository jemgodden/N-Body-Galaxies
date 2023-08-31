import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Config import *

names = []
x = []
y = []  # x, y and z coordinates of particle being plotted.
z = []
colour = []  # Colour of particle in plot.

PriGal = [[], [], []]  # List of x, y and z coordinates for position of primary galaxy during interaction.
SecGal = [[], [], []]  # List of x, y and z coordinates for position of secondary galaxy during interaction.

PriGal2 = [[], [], []]
SecGal2 = [[], [], []]

rotate = [0, 0, 0]


def option_checks():
    if centre_mid and (centre_pri or centre_sec or origin_pri or origin_sec) or \
            centre_pri and (centre_sec or origin_pri or origin_sec) or \
            centre_sec and (origin_pri or origin_sec) or \
            origin_pri and origin_sec:
        print("\nError. Please select only one centring.")
        exit(1)

    if centre_mid and not secondary_gal:
        print("\nError. Cannot centre on the middle of the interactions if there is only one galaxy.")
        exit(1)

    if centre_pri and not primary_gal:
        print("\nError. Cannot centre on the primary galaxy if there is no primary galaxy.")
        exit(1)

    if centre_sec and not secondary_gal:
        print("\nError. Cannot centre on the secondary galaxy if there is no secondary galaxy.")
        exit(1)

    if origin_pri and not primary_gal:
        print("\nError. Cannot centre the primary galaxy at the origin if there is no primary galaxy.")
        exit(1)

    if origin_sec and not secondary_gal:
        print("\nError. Cannot centre the secondary galaxy at the origin if there is no secondary galaxy.")
        exit(1)


def path_read():  # Reads information about the path of the galaxies in the interaction.
    print("\nReading files...\n\n")
    if rewind:
        file1 = open("Backwards/RWPriGalPath.txt", "r")
    else:
        file1 = open("../Forwards/PriGalPath.txt",
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
        file2 = open("../Forwards/SecGalPath.txt",
                     "r")  # Opens secondary galaxy file containing positions at all time steps.
    for line in file2:  # Goes through file line by line.
        data2 = line.strip().split()
        SecGal[0].append(float(data2[0]))
        SecGal[1].append(float(data2[1]))  # Appends each value to list to be plotted later.
        SecGal[2].append(float(data2[2]))
    file2.close()


def path_copy():
    for i in range(len(PriGal[0])):
        PriGal2[0].append(PriGal[0][i])
        PriGal2[1].append(PriGal[1][i])
        PriGal2[2].append(PriGal[2][i])  # Appends galaxy paths to editable list.
        SecGal2[0].append(SecGal[0][i])
        SecGal2[1].append(SecGal[1][i])
        SecGal2[2].append(SecGal[2][i])


def find_galaxy(galaxy_name):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def centring_mid():  # Centres the path in the images of the interaction to a central point of the two galaxies.
    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    dx1 = (x[pri] + x[sec]) / 2
    dy1 = (y[pri] + y[sec]) / 2  # Finding the mid-point of the two galaxies.
    dz1 = (z[pri] + z[sec]) / 2
    for i in range(len(x)):
        x[i] += - dx1
        y[i] += - dy1  # Altering the position of each particle in the simulation.
        z[i] += - dz1

    dx2 = []
    dy2 = []  # List of the change in position that each particle will undergo to be centred.
    dz2 = []
    for k in range(len(PriGal[0])):
        dx2.append((PriGal[0][k] + SecGal[0][k]) / 2)
        dy2.append((PriGal[1][k] + SecGal[1][k]) / 2)  # Making lists of changes in position.
        dz2.append((PriGal[2][k] + SecGal[2][k]) / 2)
    for j in range(len(PriGal[0])):
        PriGal2[0][j] += - dx2[j]
        SecGal2[0][j] += - dx2[j]
        PriGal2[1][j] += - dy2[j]
        SecGal2[1][j] += - dy2[j]  # Altering positions of each particle to centre the images.
        PriGal2[2][j] += - dz2[j]
        SecGal2[2][j] += - dz2[j]


def centring_pri():  # A function to centre the images of the interaction on the primary galaxy.
    pri = find_galaxy(pri_galaxy_name)

    dx = x[pri]
    dy = y[pri]
    dz = z[pri]
    for i in range(len(x)):
        x[i] += - dx
        y[i] += - dy  # Altering the position of each particle in the simulation.
        z[i] += - dz

    for j in range(len(PriGal[0])):
        PriGal2[0][j] += - dx
        PriGal2[1][j] += - dy
        PriGal2[2][j] += - dz  # Altering positions of each galaxy path to centre the images.
        SecGal2[0][j] += - dx
        SecGal2[1][j] += - dy
        SecGal2[2][j] += - dz


def centring_sec():  # A function to centre the images of the interaction on the secondary galaxy.
    sec = find_galaxy(sec_galaxy_name)

    dx = x[sec]
    dy = y[sec]
    dz = z[sec]
    for i in range(len(x)):
        x[i] += - dx
        y[i] += - dy  # Altering the position of each particle in the simulation.
        z[i] += - dz

    for j in range(len(PriGal[0])):
        PriGal2[0][j] += - dx
        PriGal2[1][j] += - dy
        PriGal2[2][j] += - dz  # Altering positions of each galaxy path to centre the images.
        SecGal2[0][j] += - dx
        SecGal2[1][j] += - dy
        SecGal2[2][j] += - dz


def centring_pri_origin():  # A function to keep the primary galaxy at the origin in the interaction images.
    pri = find_galaxy(pri_galaxy_name)

    dx = x[pri]
    dy = y[pri]
    dz = z[pri]
    for i in range(len(x)):
        x[i] += -dx
        y[i] += -dy
        z[i] += -dz

    for j in range(len(PriGal[0])):
        PriGal2[0][j] += - PriGal[0][j]
        SecGal2[0][j] += - PriGal[0][j]
        PriGal2[1][j] += - PriGal[1][j]
        SecGal2[1][j] += - PriGal[1][j]  # Altering positions of each particle to centre the images.
        PriGal2[2][j] += - PriGal[2][j]
        SecGal2[2][j] += - PriGal[2][j]


def centring_sec_origin():  # A function to keep the secondary galaxy at the origin in the interaction images.
    sec = find_galaxy(sec_galaxy_name)

    dx = x[sec]
    dy = y[sec]
    dz = z[sec]
    for i in range(len(x)):
        x[i] += -dx
        y[i] += -dy
        z[i] += -dz

    for j in range(len(PriGal[0])):
        PriGal2[0][j] += - SecGal[0][j]
        SecGal2[0][j] += - SecGal[0][j]
        PriGal2[1][j] += - SecGal[1][j]
        SecGal2[1][j] += - SecGal[1][j]  # Altering positions of each particle to centre the images.
        PriGal2[2][j] += - SecGal[2][j]
        SecGal2[2][j] += - SecGal[2][j]


def point_read(title):  # Reads information about all particles in the simulation.
    if rewind:
        file = open("Backwards/rimage_%.5f.txt" % title)
    else:
        file = open("Forwards/image_%.5f.txt" % title)

    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        if data[0] == pri_galaxy_name or sec_galaxy_name or pri_disk_name or sec_disk_name:
            names.append(data[0])
            x.append(float(data[2]))
            y.append(float(data[3]))  # Appends each x, y and z value to list to be plotted later.
            z.append(float(data[4]))
            colour.append(data[8])  # Appends colour of each particle to be plotted later.

    centring_choice()

    change_units()


def centring_choice():
    if centre_mid:
        centring_mid()  # For centring on the mid-point of the two galaxies.
    if centre_pri:
        centring_pri()  # For centring on the primary galaxy.
    if centre_sec:
        centring_sec()  # For centring on the secondary galaxy.
    if origin_pri:
        centring_pri_origin()  # For keeping the primary galaxy at the origin.
    if origin_sec:
        centring_sec_origin()  # For keeping the secondary galaxy at the origin.


def change_units():
    for i in range(len(x)):
        x[i] = x[i] / kpc
        y[i] = y[i] / kpc
        z[i] = z[i] / kpc

    for j in range(len(PriGal2[0])):
        PriGal2[0][j] = PriGal2[0][j] / kpc
        PriGal2[1][j] = PriGal2[1][j] / kpc
        PriGal2[2][j] = PriGal2[2][j] / kpc
        SecGal2[0][j] = SecGal2[0][j] / kpc
        SecGal2[1][j] = SecGal2[1][j] / kpc
        SecGal2[2][j] = SecGal2[2][j] / kpc


def info():  # Prints information on interaction using files made during simulation.
    total_dp1 = 0
    total_dp2 = 0
    total_part = 0
    pericentre = 0
    time_of_pericentre = 0
    if rewind:
        file = open("Backwards/RewindPericentreInfo.txt")
    else:
        file = open("../Forwards/PericentreInfo.txt")
    for line in file:  # Goes through file line by line.
        data = line.strip().split()
        total_dp1 = int(data[0])
        total_dp2 = int(data[1])
        total_part = int(data[2])
        pericentre = float(data[3])  # Pericentre of the interaction.
        time_of_pericentre = float(data[4])  # Time at which pericentre occurs during interaction.

    if primary_disk and not secondary_isolation:
        print("Primary galaxy disk particles: ", total_dp1, ".\n")
    if secondary_disk and not primary_isolation:
        print("Secondary galaxy disk particles: ", total_dp2, ".\n")
    print("Total simulation particles: ", total_part, ".\n")
    if primary_gal and secondary_gal:
        print("Pericentre: ", round(pericentre, 2), "kpc.\n")
        print("Time of pericentre: ", round(time_of_pericentre, 2), "Gyrs.\n")
    print("Number of images produced: ", (frames + 1), ".\n")
    print("Time between images: ", image_time_step, "Gyrs.\n\n")


def plot():  # Plot images of interaction.
    print("Producing images...\n")
    for i in range(frames, -1, -1):  # Plots data of each time step on a separate figure.
        title = i * image_time_step  # (frames - i)
        path_copy()
        point_read(title)
        fig = plt.figure(figsize=(6, 6))  # Creates figure.
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.set_title('t = %.2f Gyrs' % title, fontsize=10)  # Title the figure with the time of interaction.
        ax.set_xlabel('X (kpc)', fontsize=10)
        ax.set_ylabel('Y (kpc)', fontsize=10)  # Set axis labels.
        ax.set_zlabel('Z (kpc)', fontsize=10)
        ax.set_xlim(-500, 500)
        ax.set_ylim(-500, 500)  # Set axis limits.
        ax.set_zlim(-500, 500)

        if not primary_isolation and not secondary_isolation:
            if centre_pri:
                ax.plot3D(PriGal2[0], PriGal2[1], PriGal2[2], pri_path)  # Prints edited paths for when centring images on the primary galaxy.
                if secondary_gal:
                    ax.plot3D(SecGal2[0], SecGal2[1], SecGal2[2], sec_path)
            elif centre_sec:
                if primary_gal:
                    ax.plot3D(PriGal2[0], PriGal2[1], PriGal2[2], pri_path)  # Prints edited paths for when centring images on the primary galaxy.
                ax.plot3D(SecGal2[0], SecGal2[1], SecGal2[2], sec_path)
            elif origin_pri:
                ax.plot3D(PriGal2[0], PriGal2[1], PriGal2[2], pri_path)  # Plotting path of primary galaxy.
                if secondary_gal:
                    ax.plot3D(SecGal2[0], SecGal2[1], SecGal2[2], sec_path)  # Plotting path of secondary galaxy.
            elif origin_sec:
                if primary_gal:
                    ax.plot3D(PriGal2[0], PriGal2[1], PriGal2[2], pri_path)  # Plotting path of primary galaxy.
                ax.plot3D(SecGal2[0], SecGal2[1], SecGal2[2], sec_path)  # Plotting path of secondary galaxy.
            else:
                if primary_gal:
                    ax.plot3D(PriGal2[0], PriGal2[1], PriGal2[2], pri_path)  # Plotting path of primary galaxy.
                if secondary_gal:
                    ax.plot3D(SecGal2[0], SecGal2[1], SecGal2[2], sec_path)  # Plotting path of secondary galaxy.

        for j in range(len(x)):  # Plotting all objects on figure.
            ax.plot3D([x[j]], [y[j]], [z[j]], colour[j])

        PriGal2[0].clear()
        PriGal2[1].clear()
        PriGal2[2].clear()
        SecGal2[0].clear()
        SecGal2[1].clear()
        SecGal2[2].clear()
        x.clear()
        y.clear()  # Clears all the lists used to plot image, before appending new values for next image.
        z.clear()
        colour.clear()


def main():  # Calling all functions in order.

    option_checks()

    path_read()

    info()

    plot()

    plt.show()


if __name__ == '__main__':
    main()
