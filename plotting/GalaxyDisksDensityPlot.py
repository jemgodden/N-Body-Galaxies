import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Config import *

names1 = []
x1 = []
y1 = []
z1 = []
density1 = []

names2 = []
x2 = []
y2 = []
z2 = []
density2 = []

density_radius = 2.5 * kpc
final_image_time = 6.00000

pri_axis_limit = 25
sec_axis_limit = 25

point_size = 12
font_size = 18
label_size = 18

single_cb = False


def find_galaxy(names, galaxy_name):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def centring_gal(names, x, y, z, galaxy_name):
    gal = find_galaxy(names, galaxy_name)

    dx = x[gal]
    dy = y[gal]
    dz = z[gal]
    for i in range(len(x)):
        x[i] += - dx
        y[i] += - dy
        z[i] += - dz


def point_read():  # Reads information about all particles in the simulation.
    file = open("Forwards/image_%.5f.txt" % final_image_time)

    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        if data[0] == pri_galaxy_name or pri_disk_name:
            names1.append(data[0])
            x1.append(float(data[2]))
            y1.append(float(data[3]))
            z1.append(float(data[4]))
        if data[0] == sec_galaxy_name or sec_disk_name:
            names2.append(data[0])
            x2.append(float(data[2]))
            y2.append(float(data[3]))
            z2.append(float(data[4]))


def change_units(x, y, z):
    for i in range(len(x)):
        x[i] = x[i] / kpc
        y[i] = y[i] / kpc
        z[i] = z[i] / kpc


def find_separation(i, j, x, y, z):
    rx = x[i] - x[j]
    ry = y[i] - y[j]
    rz = z[i] - z[j]
    r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2))
    return r


def get_density(x, y, z, density, galaxy_name):
    gal = find_galaxy(x, galaxy_name)

    for i in range(len(x)):
        if i == gal:
            density.append(0)
        else:
            density_count = 0
            for j in range(len(x)):
                if j != gal and i != j:
                    r = find_separation(i, j, x, y, z)
                    if r < density_radius:
                        density_count += 1
            density.append(density_count)


def get_galaxy_data(names, x, y, z, density, galaxy_name):
    get_density(x, y, z, density, galaxy_name)
    centring_gal(names, x, y, z, galaxy_name)
    change_units(x, y, z)


def find_specific_separation(xi, yi, zi, xj, yj, zj):
    rx = xi - xj
    ry = yi - yj
    rz = zi - zj
    r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2))
    return r


def refine():
    pri = find_galaxy(x1, pri_galaxy_name)
    sec = find_galaxy(x2, sec_galaxy_name)

    pri_list = []
    sec_list = []

    for i in range(len(x1)):
        if i == pri:
            continue
        else:
            pri_r = find_specific_separation(x1[pri], y1[pri], z1[pri], x1[i], y1[i], z1[i])
            sec_r = find_specific_separation(x2[sec], y2[sec], z2[sec], x1[i], y1[i], z1[i])
            if pri_r > 30 and sec_r > 30:
                pri_list.append(i)

    for j in range(len(x2)):
        if j == sec:
            continue
        else:
            pri_r = find_specific_separation(x1[pri], y1[pri], z1[pri], x2[j], y2[j], z2[j])
            sec_r = find_specific_separation(x2[sec], y2[sec], z2[sec], x2[j], y2[j], z2[j])
            if pri_r > 30 and sec_r > 30:
                sec_list.append(j)

    for m in range(len(pri_list)-1, 0, -1):
        names1.pop(pri_list[m])
        x1.pop(pri_list[m])
        y1.pop(pri_list[m])
        z1.pop(pri_list[m])
        density1.pop(pri_list[m])

    for n in range(len(sec_list)-1, 0, -1):
        names2.pop(pri_list[n])
        x2.pop(pri_list[n])
        y2.pop(pri_list[n])
        z2.pop(pri_list[n])
        density2.pop(pri_list[n])


def plot():  # Plot images of interaction.
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

    if primary_gal:
        ax1.set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        ax1.set_ylabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        ax1.set_xlim(-pri_axis_limit, pri_axis_limit)
        ax1.set_ylim(-pri_axis_limit, pri_axis_limit)
        ax1.tick_params(labelsize=label_size)

        ax1.scatter(x1, y1, c=density1, cmap='plasma', s=point_size)
        ax2.set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        ax2.set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        ax2.set_xlim(-pri_axis_limit, pri_axis_limit)
        ax2.set_ylim(-pri_axis_limit, pri_axis_limit)
        ax2.tick_params(labelsize=label_size)

        ax2.scatter(x1, z1, c=density1, cmap='plasma', s=point_size)
        ax3.set_xlabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        ax3.set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        ax3.set_xlim(-pri_axis_limit, pri_axis_limit)
        ax3.set_ylim(-pri_axis_limit, pri_axis_limit)
        ax3.tick_params(labelsize=label_size)
        im3 = ax3.scatter(y1, z1, c=density1, cmap='plasma', s=point_size)
        divider = make_axes_locatable(ax3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        cb1 = fig.colorbar(im3, cax=cax1)
        cb1.set_label(label='Nearest Neighbour Density', size=16)
        cb1.ax.tick_params(labelsize=16)

    if secondary_gal:
        ax4.set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        ax4.set_ylabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        ax4.set_xlim(-sec_axis_limit, sec_axis_limit)
        ax4.set_ylim(-sec_axis_limit, sec_axis_limit)
        ax4.tick_params(labelsize=label_size)

        ax4.scatter(x2, y2, c=density2, cmap='plasma', s=point_size)
        ax5.set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        ax5.set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        ax5.set_xlim(-sec_axis_limit, sec_axis_limit)
        ax5.set_ylim(-sec_axis_limit, sec_axis_limit)
        ax5.tick_params(labelsize=label_size)

        ax5.scatter(x2, z2, c=density2, cmap='plasma', s=point_size)
        ax6.set_xlabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        ax6.set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        ax6.set_xlim(-sec_axis_limit, sec_axis_limit)
        ax6.set_ylim(-sec_axis_limit, sec_axis_limit)
        ax6.tick_params(labelsize=label_size)
        im6 = ax6.scatter(y2, z2, c=density2, cmap='plasma', s=point_size)
        divider = make_axes_locatable(ax6)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        cb2 = fig.colorbar(im6, cax=cax1)
        cb2.set_label(label='Nearest Neighbour Density', size=16)
        cb2.ax.tick_params(labelsize=16)

    fig.tight_layout(pad=-2.0)
    plt.show()


def single_cb_plot():

    fig, axs = plt.subplots(2, 3)

    if primary_gal:
        axs[0, 0].set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 0].set_ylabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 0].set_xlim(-pri_axis_limit, pri_axis_limit)
        axs[0, 0].set_ylim(-pri_axis_limit, pri_axis_limit)
        axs[0, 0].tick_params(labelsize=label_size)
        axs[0, 0].scatter(x1, y1, c=density1, cmap='plasma', s=point_size)

        axs[0, 1].set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 1].set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 1].set_xlim(-pri_axis_limit, pri_axis_limit)
        axs[0, 1].set_ylim(-pri_axis_limit, pri_axis_limit)
        axs[0, 1].tick_params(labelsize=label_size)
        axs[0, 1].scatter(x1, z1, c=density1, cmap='plasma', s=point_size)

        axs[0, 2].set_xlabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 2].set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        axs[0, 2].set_xlim(-pri_axis_limit, pri_axis_limit)
        axs[0, 2].set_ylim(-pri_axis_limit, pri_axis_limit)
        axs[0, 2].tick_params(labelsize=label_size)
        cbp = axs[0, 2].scatter(y1, z1, c=density1, cmap='plasma', s=point_size)

    if secondary_gal:
        axs[1, 0].set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 0].set_ylabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 0].set_xlim(-sec_axis_limit, sec_axis_limit)
        axs[1, 0].set_ylim(-sec_axis_limit, sec_axis_limit)
        axs[1, 0].tick_params(labelsize=label_size)
        axs[1, 0].scatter(x2, y2, c=density2, cmap='plasma', s=point_size)

        axs[1, 1].set_xlabel('X $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 1].set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 1].set_xlim(-sec_axis_limit, sec_axis_limit)
        axs[1, 1].set_ylim(-sec_axis_limit, sec_axis_limit)
        axs[1, 1].tick_params(labelsize=label_size)
        axs[1, 1].scatter(x2, z2, c=density2, cmap='plasma', s=point_size)

        axs[1, 2].set_xlabel('Y $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 2].set_ylabel('Z $(kpc)$', fontsize=font_size, weight='bold')
        axs[1, 2].set_xlim(-sec_axis_limit, sec_axis_limit)
        axs[1, 2].set_ylim(-sec_axis_limit, sec_axis_limit)
        axs[1, 2].tick_params(labelsize=label_size)
        axs[1, 2].scatter(y2, z2, c=density2, cmap='plasma', s=point_size)

    plt.subplots_adjust(wspace=0.5, hspace=0.3)

    cb = fig.colorbar(cbp, ax=axs)
    cb.set_label(label='Nearest Neighbour Density', size=18)
    cb.ax.tick_params(labelsize=18)

    plt.show()


def main():  # Calling all functions in order.
    point_read()

    if primary_gal:
        get_galaxy_data(names1, x1, y1, z1, density1, pri_galaxy_name)
    if secondary_gal:
        get_galaxy_data(names2, x2, y2, z2, density2, sec_galaxy_name)

    refine()

    if single_cb:
        single_cb_plot()
    else:
        plot()


if __name__ == '__main__':
    main()
