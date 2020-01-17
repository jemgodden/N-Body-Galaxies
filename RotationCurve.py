import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

x = []
y = []  # x, y and z coordinates of particles in disk.
z = []

vx = []
vy = []  # x, y and z velocities of particles in disk.
vz = []

r = []  # List of distances from centre of galaxy to particle.
v = []  # List of velocities of particles.


def path_read():  # Reads information about the path of the galaxies in the interaction.
    print("\nReading files...\n\n")
    if primary_isolation:
        file1 = open("Primary_Galaxy.txt", "r")
        for line in file1:  # Goes through file line by line.
            data1 = line.strip().split()
            x.append(float(data1[2]))
            y.append(float(data1[3]))  # Appends each x, y and z position to list.
            z.append(float(data1[4]))
            vx.append(float(data1[5]))
            vy.append(float(data1[6]))  # Appends each x, y and z velocity to list.
            vz.append(float(data1[7]))
        file1.close()

    if secondary_isolation:
        file2 = open("Secondary_Galaxy.txt", "r")
        for line in file2:  # Goes through file line by line.
            data2 = line.strip().split()
            x.append(float(data2[2]))
            y.append(float(data2[3]))  # Appends each x, y and z position to list.
            z.append(float(data2[4]))
            vx.append(float(data2[5]))
            vy.append(float(data2[6]))  # Appends each x, y and z velocity to list.
            vz.append(float(data2[7]))
        file2.close()


def calc_data():
    for i in range(len(x)-1):
        d = math.sqrt((x[0]-x[i+1])**2 + (y[0]-y[i+1])**2 + (z[0]-z[i+1])**2)  # Finds distance from centre of galaxy to particle.
        r.append(d/kpc)
        u = math.sqrt((vx[i+1])**2 + (vy[i+1])**2 + (vz[i+1])**2)    # Finds velocity of particle.
        v.append(u/km_s)


def plot():
    print("Producing images...\n")
    plt.plot(r, v, 'k-')  # Plots distcance from centre against particle velocity.
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("Distance from centre of galaxy (kpc)")
    plt.ylabel("Radial velocity (km/s)")


def main():  # Calling all functions in order.

    if not primary_isolation and not secondary_isolation or primary_isolation and secondary_isolation:
        print("\nError, can only plot a galaxy rotation curve when there is a galaxy simulated in isolation.\n")
        exit(1)

    path_read()

    calc_data()

    plot()

    plt.show()


if __name__ == '__main__':
    main()
