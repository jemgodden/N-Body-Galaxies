import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from IntData import *

pri_names = []
pri_x = []
pri_y = []
pri_z = []
sec_names = []
sec_x = []
sec_y = []
sec_z = []

pri_rad = []
sec_rad = []

initial_image_time = 0.00000
final_image_time = 1.50000

pri_axis_limit = 25
sec_axis_limit = 25

point_size = 20
font_size = 18
label_size = 18

x_los = [0, 4.113]
y_los = [0, -8.853]
z_los = [0, -2.178]


def find_galaxy(galaxy_name, names):  # Finds position of a galaxy in the list of bodies.
    position = 0
    for i in range(len(names)):
        if names[i] == galaxy_name:
            position = i
            break
    return position


def centring_particles():
    pri = find_galaxy(pri_galaxy_name, pri_names)
    sec = find_galaxy(sec_galaxy_name, sec_names)

    dx1 = (pri_x[pri] + sec_x[sec]) / 2
    dy1 = (pri_y[pri] + sec_y[sec]) / 2
    dz1 = (pri_z[pri] + sec_z[sec]) / 2
    for i in range(len(pri_x)):
        pri_x[i] += - dx1
        sec_x[i] += - dx1
        pri_y[i] += - dy1
        sec_y[i] += - dy1
        pri_z[i] += - dz1
        sec_z[i] += - dz1


def point_read(image_time):  # Reads information about all particles in the simulation.
    file = open("Forwards/image_%.5f.txt" % image_time)

    for line in file:  # Reads the file line by line.
        data = line.strip().split()
        if data[0] == pri_galaxy_name or data[0] == pri_disk_name:
            pri_names.append(data[0])
            pri_x.append(float(data[2]))
            pri_y.append(float(data[3]))
            pri_z.append(float(data[4]))
        elif data[0] == sec_galaxy_name or sec_disk_name:
            sec_names.append(data[0])
            sec_x.append(float(data[2]))
            sec_y.append(float(data[3]))
            sec_z.append(float(data[4]))


def change_units():
    for i in range(len(pri_x)):
        pri_x[i] = pri_x[i] / kpc
        sec_x[i] = sec_x[i] / kpc
        pri_y[i] = pri_y[i] / kpc
        sec_y[i] = sec_y[i] / kpc
        pri_z[i] = pri_z[i] / kpc
        sec_z[i] = sec_z[i] / kpc


def find_separation(i, j, x, y, z):
    rx = x[i] - x[j]
    ry = y[i] - y[j]
    rz = z[i] - z[j]
    r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2)) / kpc
    return r


def find_disk_particle_radii(galaxy_name, names, x, y, z, rad):
    gal = find_galaxy(galaxy_name, names)
    for i in range(len(x)):
        if i == gal:
            rad.append(0)
        else:
            rad.append(find_separation(gal, i, x, y, z))


def get_galaxy_data():
    centring_particles()
    change_units()


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    pri_plot = ax.scatter(pri_x, pri_y, pri_z, c=pri_rad, cmap='Blues', s=point_size, depthshade=False)
    sec_plot = ax.scatter(sec_x, sec_y, sec_z, c=sec_rad, cmap='Reds', s=point_size, depthshade=False)
    # ax.plot(x_los, y_los, z_los, 'k-')
    ax.set_axis_off()

    pri_cb = fig.colorbar(pri_plot)
    pri_cb.set_label(label='Initial Radius from Centre of NGC5257 $(kpc)$', size=18)
    pri_cb.ax.tick_params(labelsize=18)
    sec_cb = fig.colorbar(sec_plot)
    sec_cb.set_label(label='Initial radius from centre of NGC5258 $(kpc)$', size=18)
    sec_cb.ax.tick_params(labelsize=18)
    plt.show()


def main():  # Calling all functions in order.
    point_read(initial_image_time)

    find_disk_particle_radii(pri_galaxy_name, pri_names, pri_x, pri_y, pri_z, pri_rad)
    find_disk_particle_radii(sec_galaxy_name, sec_names, sec_x, sec_y, sec_z, sec_rad)

    pri_names.clear()
    pri_x.clear()
    pri_y.clear()
    pri_z.clear()
    sec_names.clear()
    sec_x.clear()
    sec_y.clear()
    sec_z.clear()

    point_read(final_image_time)

    get_galaxy_data()

    plot()


if __name__ == '__main__':
    main()
