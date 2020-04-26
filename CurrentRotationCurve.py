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


def get_vector_magnitude(v):
    return math.sqrt((v[0] ** 2) + (v[1] ** 2) + (v[2] ** 2))


def dot_product(v1, v2):
    x = 0
    for i in range(len(v1)):
        x += v1[i] * v2[i]
    return x


def angle_between_vectors(v1, v2):
    v1_mag = get_vector_magnitude(v1)
    v2_mag = get_vector_magnitude(v2)
    dot_prod = dot_product(v1, v2)
    return math.acos(dot_prod / (v1_mag * v2_mag))


def cross_product(v1, v2):
    n = []
    n.append((v2[2] * v1[1]) - (v2[1] * v1[2]))
    n.append((v2[0] * v1[2]) - (v2[2] * v1[0]))
    n.append((v2[1] * v1[0]) - (v2[0] * v1[1]))
    return n


def get_radial_velocity(uxyz, dxyz, n):
    x_prod = cross_product(dxyz, n)
    theta = angle_between_vectors(uxyz, x_prod)
    vxyz = []
    for i in range(len(uxyz)):
        vxyz.append(uxyz[i] * math.cos(theta))
    return get_vector_magnitude(vxyz)


def find_radial_velocity(galaxy_name, names, x, y, z, vx, vy, vz, r, v):
    if galaxy_name == pri_galaxy_name:
        max_r = dr1
        n = norm_spin1
    elif galaxy_name == sec_galaxy_name:
        max_r = dr2
        n = norm_spin2

    gal = find_galaxy(galaxy_name, names)

    for i in range(len(x)-1):
        if i == gal:
            continue
        else:
            dxyz = [x[i] - x[gal], y[i] - y[gal], z[i] - z[gal]]
            d = math.sqrt(dxyz[0] ** 2 + dxyz[1] ** 2 + dxyz[2] ** 2)
            if d < max_r:
                r.append(d/kpc)

                uxyz = [vx[i] - vx[gal], vy[i] - vy[gal], vz[i] - vz[gal]]
                u = get_radial_velocity(uxyz, dxyz, n)
                v.append(u/km_s)
            else:
                continue


def func(x_plot, a, b, c):
    f = []
    for i in range(len(x_plot)):
        f.append(a * np.exp(-b * x_plot[i]) + c)
    return f


def plot_rot_curve(galaxy_id, r, v):
    r, v = (list(t) for t in zip(*sorted(zip(r, v))))

    # bins = 15
    # if galaxy_id == primary:
    #     d = dr1 / bins
    # elif galaxy_id == secondary:
    #     d = dr2 / bins
    # v_avg = []
    # r_bin = []
    # for j in range(bins):
    #     n = 0
    #     v_tot = 0
    #     for i in range(len(r)):
    #         if (j * d / kpc) < r[i] <= ((j + 1) * d / kpc):
    #             n += 1
    #             v_tot += v[i]
    #         else:
    #             continue
    #     if n == 0:
    #         continue
    #     else:
    #         v_avg.append(v_tot / n)
    #         r_bin.append((j + 0.5) * d / kpc)

    popt1, pcov1 = curve_fit(func, r, v, p0=(1, 1, 1))
    # popt2, pcov2 = curve_fit(func, r_bin, v_avg, p0=(1, 1, 1))

    plt.figure(figsize=(6, 6))
    plt.plot(r, v, 'b.', markersize=18)
    plt.plot(r, func(r, *popt1), 'r-', linewidth=8)
    # plt.plot(r_bin, v_avg, 'r-', linewidth=8)
    # plt.plot(r_bin, func(r_bin, *popt2), 'r-', linewidth=8)
    plt.gca().set_ylim(bottom=0)
    if galaxy_id == primary:
        plt.xlabel("Distance from Centre of NGC 5257 $(kpc)$", weight='bold', fontsize=28)
    elif galaxy_id == secondary:
        plt.xlabel("Distance from Centre of NGC 5258 $(kpc)$", weight='bold', fontsize=28)
    else:
        plt.xlabel("Distance from centre of galaxy $(kpc)$", weight='bold', fontsize=28)
    plt.ylabel("Radial velocity $(km/s)$", weight='bold', fontsize=28)
    plt.tick_params(labelsize=26)

    plt.show()


def main():  # Calling all functions in order.

    file_read()

    find_radial_velocity(pri_galaxy_name, pri_names, pri_x, pri_y, pri_z, pri_vx, pri_vy, pri_vz, pri_r, pri_v)
    find_radial_velocity(sec_galaxy_name, sec_names, sec_x, sec_y, sec_z, sec_vx, sec_vy, sec_vz, sec_r, sec_v)

    plot_rot_curve(primary, pri_r, pri_v)
    plot_rot_curve(secondary, sec_r, sec_v)


if __name__ == '__main__':
    main()
