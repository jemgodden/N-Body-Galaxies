import math
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from Config import *

names = []
x = []
y = []  # x, y and z coordinates of particles in disk.
z = []

vx = []
vy = []  # x, y and z velocities of particles in disk.
vz = []

r = []  # List of distances from centre of galaxy to particle.
v = []  # List of velocities of particles.


def file_read():  # Reads information about the path of the galaxies in the interaction.
    print("\nReading files...\n\n")
    if primary_isolation:
        file1 = open("../../data/Primary_Galaxy.txt", "r")
        for line in file1:  # Goes through file line by line.
            data1 = line.strip().split()
            names.append(data1[0])
            x.append(float(data1[2]))
            y.append(float(data1[3]))  # Appends each x, y and z position to list.
            z.append(float(data1[4]))
            vx.append(float(data1[5]))
            vy.append(float(data1[6]))  # Appends each x, y and z velocity to list.
            vz.append(float(data1[7]))
        file1.close()

    if secondary_isolation:
        file2 = open("../../data/Secondary_Galaxy.txt", "r")
        for line in file2:  # Goes through file line by line.
            data2 = line.strip().split()
            names.append(data2[0])
            x.append(float(data2[2]))
            y.append(float(data2[3]))  # Appends each x, y and z position to list.
            z.append(float(data2[4]))
            vx.append(float(data2[5]))
            vy.append(float(data2[6]))  # Appends each x, y and z velocity to list.
            vz.append(float(data2[7]))
        file2.close()


def find_galaxy(galaxy_name):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def calc_data():
    gal = 0
    if primary_isolation:
        gal = find_galaxy(pri_galaxy_name)
    if secondary_isolation:
        gal = find_galaxy(sec_galaxy_name)

    for i in range(len(x)-1):
        if gal == i:
            continue
        else:
            d = math.sqrt((x[gal]-x[i])**2 + (y[gal]-y[i])**2 + (z[gal]-z[i])**2)
            r.append(d/kpc)
            u = math.sqrt((vx[i])**2 + (vy[i])**2 + (vz[i])**2)    # Finds velocity of particle.
            v.append(u/km_s)


def plot_rot_curve():
    plt.figure(figsize=(6, 6))
    plt.plot(r, v, 'k.', markersize=18)  # Plots distance from centre against particle velocity.
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("Distance from Centre of NGC5258 $(kpc)$", weight='bold', fontsize=28)
    plt.ylabel("Radial velocity $(km/s)$", weight='bold', fontsize=28)
    plt.tick_params(labelsize=26)

    plt.show()


def func(x_plot, a, b, c):
    f = []
    for i in range(len(x_plot)):
        f.append(a * math.exp(-b * x_plot[i]) + c)
    return f


def plot_distributions():
    bins = 50
    extra = (0.05 * bins) // 1
    d = dr1 / bins
    dist = []
    n = [0] * bins
    for i in range(len(r)):
        for j in range(bins):
            if (j + extra) * d < r[i] * kpc <= (j + extra + 1) * d:
                n[j] += 1
            else:
                continue

    for k in range(len(n)):
        dist.append(((k + extra + 0.5) * d) / kpc)
        n[k] = n[k]

    coefs = poly.polyfit(dist, n, 4)
    ffit = poly.polyval(dist, coefs)

    plt.figure(figsize=(6, 6))
    plt.plot(dist, n, 'k.')
    plt.plot(dist, ffit, 'k-')
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("Distance from centre of galaxy (kpc)")
    plt.ylabel("Number of particles")

    plt.show()

    for q in range(len(n)):
        area = math.pi * ((d / kpc) ** 2) * (((q + 1) ** 2) - (q ** 2))
        n[q] = n[q] / area

    popt, pcov = curve_fit(func, dist, n)

    plt.figure(figsize=(6, 6))
    plt.plot(dist, n, 'k.')
    plt.plot(dist, func(dist, *popt), 'k-')
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("Distance from centre of galaxy (kpc)")
    plt.ylabel("Surface density of particles (particles $kpc^{-2}$)")

    plt.show()


def main():  # Calling all functions in order.

    if not primary_isolation and not secondary_isolation or primary_isolation and secondary_isolation:
        print("\nError, can only plot a galaxy rotation curve when there is a galaxy simulated in isolation.\n")
        exit(1)

    file_read()

    calc_data()

    plot_rot_curve()

    plot_distributions()

    # plt.show()


if __name__ == '__main__':
    main()
