import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Config import *

names = []
x = []
y = []
z = []
density = []

density_radius = 2.5 * kpc
final_image_time = 1.50000

point_size = 20
font_size = 18
label_size = 18

x_los = [0, 41.13]
y_los = [0, -88.53]
z_los = [0, -21.78]


def find_galaxy(galaxy_name):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def centring_particles():
    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    dx1 = (x[pri] + x[sec]) / 2
    dy1 = (y[pri] + y[sec]) / 2
    dz1 = (z[pri] + z[sec]) / 2
    for i in range(len(x)):
        x[i] += - dx1
        y[i] += - dy1
        z[i] += - dz1


def point_read():  # Reads information about all particles in the simulation.
    file = open("Forwards/image_%.5f.txt" % final_image_time)

    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        if data[0] == pri_galaxy_name or sec_galaxy_name or pri_disk_name or sec_disk_name:
            names.append(data[0])
            x.append(float(data[2]))
            y.append(float(data[3]))
            z.append(float(data[4]))


def change_units():
    for i in range(len(x)):
        x[i] = x[i] / kpc
        y[i] = y[i] / kpc
        z[i] = z[i] / kpc


def find_separation(i, j):
    rx = x[i] - x[j]
    ry = y[i] - y[j]
    rz = z[i] - z[j]
    r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2))
    return r


def get_density():
    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    for i in range(len(x)):
        if i == pri or i == sec:
            density.append(0)
        else:
            density_count = 0
            for j in range(len(x)):
                if j != pri and j != sec and i != j:
                    r = find_separation(i, j)
                    if r < density_radius:
                        density_count += 1
            density.append(density_count)


def get_galaxy_data():
    get_density()
    centring_particles()
    change_units()


def refine():
    pri = find_galaxy(pri_galaxy_name)
    sec = find_galaxy(sec_galaxy_name)

    pop_list = []

    for i in range(len(x)):
        if i == pri or i == sec:
            continue
        else:
            pri_r = find_separation(pri, i)
            sec_r = find_separation(sec, i)
            if pri_r > 300 and sec_r > 300:
                pop_list.append(i)

    for m in range(len(pop_list)-1, 0, -1):
        names.pop(pop_list[m])
        x.pop(pop_list[m])
        y.pop(pop_list[m])
        z.pop(pop_list[m])
        density.pop(pop_list[m])


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    p = ax.scatter(x, y, z, c=density, cmap='plasma', s=point_size)
    ax.set_axis_off()
    # ax.plot(x_los, y_los, z_los, 'k-')  # Points towards Earth to help line up line of sight.
    cb = fig.colorbar(p)
    cb.set_label(label='Nearest Neighbour Density', size=18)
    cb.ax.tick_params(labelsize=18)
    plt.show()


def main():  # Calling all functions in order.
    point_read()

    get_galaxy_data()

    refine()

    plot()


if __name__ == '__main__':
    main()
