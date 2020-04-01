import math
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from IntData import *

pri_names = []
pri_x = []
pri_y = []
pri_z = []
pri_vx = []
pri_vy = []
pri_vz = []

sec_names = []
sec_x = []
sec_y = []
sec_z = []
sec_vx = []
sec_vy = []
sec_vz = []

pri_r = []
pri_v = []

sec_r = []
sec_v = []


def file_read():  # Reads information about the path of the galaxies in the interaction.
    file = open("Forwards/image_1.50000.txt", "r")
    for line in file:  # Goes through file line by line.
        data = line.strip().split()
        if data[0] == pri_galaxy_name or pri_disk_name:
            pri_names.append(data[0])
            pri_x.append(float(data[2]))
            pri_y.append(float(data[3]))
            pri_z.append(float(data[4]))
            pri_vx.append(float(data[5]))
            pri_vy.append(float(data[6]))
            pri_vz.append(float(data[7]))
        if data[0] == sec_galaxy_name or sec_disk_name:
            sec_names.append(data[0])
            sec_x.append(float(data[2]))
            sec_y.append(float(data[3]))
            sec_z.append(float(data[4]))
            sec_vx.append(float(data[5]))
            sec_vy.append(float(data[6]))
            sec_vz.append(float(data[7]))
    file.close()


def find_galaxy(galaxy_name, names):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def calc_data(galaxy_name, names, x, y, z, vx, vy, vz, r, v):
    gal = find_galaxy(galaxy_name, names)

    for i in range(len(x)-1):
        if i == gal:
            continue
        else:
            d = math.sqrt((x[gal]-x[i])**2 + (y[gal]-y[i])**2 + (z[gal]-z[i])**2)
            if d < (30 * kpc):
                r.append(d/kpc)
                u = math.sqrt((vx[i])**2 + (vy[i])**2 + (vz[i])**2)
                v.append(u/km_s)
            else:
                continue


def plot_rot_curve(galaxy_id, r, v):
    plt.figure(figsize=(6, 6))
    plt.plot(r, v, 'k.', markersize=18)  # Plots distance from centre against particle velocity.
    plt.gca().set_ylim(bottom=0)
    if galaxy_id == primary:
        plt.xlabel("Distance from centre of primary galaxy $(kpc)$", weight='bold', fontsize=28)
    elif galaxy_id == secondary:
        plt.xlabel("Distance from centre of secondary galaxy $(kpc)$", weight='bold', fontsize=28)
    else:
        plt.xlabel("Distance from centre of galaxy $(kpc)$", weight='bold', fontsize=28)
    plt.ylabel("Radial velocity $(km/s)$", weight='bold', fontsize=28)
    plt.tick_params(labelsize=26)

    plt.show()


def main():  # Calling all functions in order.

    file_read()

    calc_data(pri_galaxy_name, pri_names, pri_x, pri_y, pri_z, pri_vx, pri_vy, pri_vz, pri_r, pri_v)
    calc_data(sec_galaxy_name, sec_names, sec_x, sec_y, sec_z, sec_vx, sec_vy, sec_vz, sec_r, sec_v)

    plot_rot_curve(primary, pri_r, pri_v)
    plot_rot_curve(secondary, sec_r, sec_v)


if __name__ == '__main__':
    main()
