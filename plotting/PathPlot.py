import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Config import *

PriPath = [[], [], []]
SecPath = [[], [], []]

Centre = [[0, 0], [0, 0], [0, 0]]

fraction = 6 / 6

SecGalCentre = True
TwoD = True


def path_read():
    max = round(30000 * fraction)
    i = 0
    j = 0

    file1 = open("../Forwards/PriGalPath.txt", "r")
    # file1 = open("Backwards/RWPriGalPath.txt", "r")
    for line in file1:
        data1 = line.strip().split()
        PriPath[0].append(float(data1[0]))
        PriPath[1].append(float(data1[1]))
        PriPath[2].append(float(data1[2]))
        i += 1
        if i > max:
            break
    file1.close()

    file2 = open("../Forwards/SecGalPath.txt", "r")
    # file2 = open("Backwards/RWSecGalPath.txt", "r")
    for line in file2:
        data2 = line.strip().split()
        SecPath[0].append(float(data2[0]))
        SecPath[1].append(float(data2[1]))
        SecPath[2].append(float(data2[2]))
        j += 1
        if j > max:
            break
    file2.close()


def find_vectors(v1, v2):
    for i in range(len(PriPath)):
        v1.append(PriPath[i][-1] - PriPath[i][0])
    for i in range(len(PriPath)):
        v2.append(SecPath[i][-1] - PriPath[i][0])


def find_plane_normal(v1, v2, n):
    n.append((v1[1]*v2[2])-(v1[2]*v2[1]))
    n.append(-((v1[0] * v2[2]) - (v1[2] * v2[0])))
    n.append((v1[0] * v2[1]) - (v1[1] * v2[0]))


def move_particle(alpha, beta, x, y, z):
    xyz_temp = [x, y, z]

    x = xyz_temp[0] * math.cos(beta) + (xyz_temp[1] * math.cos(alpha) + xyz_temp[2] * math.sin(alpha)) * math.sin(beta)
    y = - xyz_temp[0] * math.sin(beta) + (xyz_temp[1] * math.cos(alpha) + xyz_temp[2] * math.sin(alpha)) * math.cos(beta)
    z = - xyz_temp[1] * math.sin(alpha) + xyz_temp[2] * math.cos(alpha)


def transform_plane(normal):
    if normal[2] == -1:
        alpha = math.pi
    else:
        alpha = 2 * math.atan((math.sqrt((normal[1] ** 2) + (normal[2] ** 2) - 1) - normal[1]) / (normal[2] + 1))
    beta = math.atan(- normal[0] / ((normal[1] * math.cos(alpha)) + (normal[2] * math.sin(alpha))))

    for i in range(len(PriPath[0])):
        move_particle(alpha, beta, PriPath[0][i], PriPath[1][i], PriPath[2][i])

    for j in range(len(SecPath[0])):
        move_particle(alpha, beta, SecPath[0][j], SecPath[1][j], SecPath[2][j])


def make_2D():
    vector1 = []
    vector2 = []
    normal = []

    find_vectors(vector1, vector2)
    find_plane_normal(vector1, vector2, normal)
    transform_plane(normal)


def centring_sec_origin():
    for i in range(len(PriPath[0])):
        PriPath[0][i] += - SecPath[0][i]
        SecPath[0][i] += - SecPath[0][i]
        PriPath[1][i] += - SecPath[1][i]
        SecPath[1][i] += - SecPath[1][i]
        PriPath[2][i] += - SecPath[2][i]
        SecPath[2][i] += - SecPath[2][i]


def change_units():
    for i in range(len(PriPath)):
        for j in range(len(PriPath[0])):
            PriPath[i][j] = PriPath[i][j] / kpc
            SecPath[i][j] = SecPath[i][j] / kpc


def no_axis_3D_plot():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_axis_off()
    if SecGalCentre:
        ax.plot(Centre[0], Centre[1], Centre[2], 'r.')
    ax.plot(PriPath[0], PriPath[1], PriPath[2], pri_path)
    ax.plot(SecPath[0], SecPath[1], SecPath[2], sec_path)
    plt.show()


def plot_2D():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    set1, = ax.plot(PriPath[0], PriPath[1], color='blue', linewidth=10, label='NGC 5257')
    if SecGalCentre:
        set2 = ax.scatter(0, 0, color='red', s=250, label='NGC 5258')
        lgnd = ax.legend(handles=[set1, set2], fontsize=24, markerscale=3)
        lgnd.legendHandles[1]._sizes = [220]
    else:
        ax.plot(SecPath[0], SecPath[1], color='red', linewidth=10, label='NGC 5258')
        ax.legend(fontsize=28, markerscale=3)
    ax.set_xlabel('X - In Plane of Interaction $(kpc)$', fontsize=30, weight='bold')
    ax.set_ylabel('Y - In Plane of Interaction $(kpc)$', fontsize=30, weight='bold')
    ax.tick_params(labelsize=28)

    plt.show()


def main():
    path_read()

    if TwoD:
        make_2D()
        if SecGalCentre:
            centring_sec_origin()
        change_units()
        plot_2D()
    else:
        if SecGalCentre:
            centring_sec_origin()
        no_axis_3D_plot()


if __name__ == '__main__':
    main()
