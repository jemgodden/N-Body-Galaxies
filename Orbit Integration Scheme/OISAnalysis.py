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

# pericentre = []
# time_of_pericentre = []

pericentre1 = []
time_of_pericentre1 = []
pericentre2 = []
time_of_pericentre2 = []
pericentre3 = []
time_of_pericentre3 = []
# pericentre4 = []
# time_of_pericentre4 = []

separation = []
ratio = []
usable_lines = []
input_file_name = "Automation_Data7.txt"
output_file_name = "Adequate_Data7.txt"


def read_file(file_name, pericentre, time_of_pericentre):
    file = open(file_name, "r")
    line_no = 0
    for line in file:
        line_no += 1
        if line[0] == '#':
            continue
        else:
            data = line.strip().split("\t")
            if float(data[12]) < 0.1 or -float(data[13]) > 1.4:
                continue
            else:
                pericentre.append(float(data[12]))
                time_of_pericentre.append(-float(data[13]))
                separation.append(float(data[14]))
                # ratio.append(float(data[18]))
                # if 0.2 < -float(data[13]) < 0.3 and float(data[12]) < 50:
                #     usable_lines.append(line_no)
    file.close()


def plot_pericentre_info():
    fig, ax = plt.subplots()
    # ax.plot(pericentre, time_of_pericentre, 'k.', label="Simulation Set")
    set1, = ax.plot(pericentre1, time_of_pericentre1, 'b.', label="Simulation Set 1")
    set2, = ax.plot(pericentre2, time_of_pericentre2, 'c.', label="Simulation Set 2")
    set3, = ax.plot(pericentre3, time_of_pericentre3, 'g.', label="Simulation Set 3")
    set4 = ax.scatter(37.73, 0.2234, c='y', s=100, label="Result Used", zorder=2.5)
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    set0 = ax.fill_between([10, 40], [0.3, 0.3], [0.2, 0.2], fc='red', alpha=0.25, label="Reasonable Results", zorder=3)
    ax.set_xlabel("Pericentre Separation $(kpc)$", fontsize=25, weight='bold')
    ax.set_ylabel("Lookback Time $(Gyrs)$", fontsize=25, weight='bold')
    lgnd = ax.legend(handles=[set0, set1, set2, set3, set4], fontsize=24, markerscale=3, loc='lower left')
    lgnd.legendHandles[4]._sizes = [80]
    ax.tick_params(labelsize=22)


def plot_separation_pericentre_info():
    plt.figure(figsize=(6, 6))
    plt.plot(ratio, pericentre, 'k.')
    plt.gca().set_ylim([0, None])
    plt.xlabel("Ratio")
    plt.ylabel("Pericentre (kpc)")


def plot_separation_time_info():
    plt.figure(figsize=(6, 6))
    plt.plot(ratio, time_of_pericentre, 'k.')
    plt.gca().set_ylim([0, None])
    plt.xlabel("Ratio")
    plt.ylabel("Lookback Time of Pericentre (Gyrs)")


def plot_pericentre_info_cb():
    fig, ax = plt.subplots()
    p = ax.scatter(pericentre, time_of_pericentre, c=ratio, cmap='plasma', s=8)
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    ax.fill_between([10, 40], [0.3, 0.3], [0.2, 0.2], fc='red', alpha=0.5)
    ax.set_xlabel("Pericentre $(kpc)$", fontsize=25, weight='bold')
    ax.set_ylabel("Lookback Time $(Gyrs)$", fontsize=25, weight='bold')
    ax.tick_params(labelsize=22)
    cb = fig.colorbar(p)
    cb.set_label(label='Ratio', size=18)
    cb.ax.tick_params(labelsize=18)


def round_up_time(value):
    value = 100 * value
    new_value = math.ceil(value)

    return new_value / 100


def find_min_pericentre_time():
    input_file = open(input_file_name, "r")
    output_file = open(output_file_name, "w+")
    for line in input_file:
        if line[0] == '#':
            continue
        else:
            data = line.strip().split("\t")
            if (37.73 < float(data[12]) < 37.734) and (0.22338 < -float(data[13]) < 0.22342):
                print(line)
                output_file.write(line)
    input_file.close()
    output_file.close()


def main():  # Calling all functions in order.

    # read_file("Automation_Data7.txt", pericentre, time_of_pericentre)

    read_file("Automation_Data5.txt", pericentre1, time_of_pericentre1)
    read_file("Automation_Data6.txt", pericentre2, time_of_pericentre2)
    read_file("Automation_Data7.txt", pericentre3, time_of_pericentre3)
    # read_file("Initial_Conditions_Automation/Automation_Data4.txt", pericentre4, time_of_pericentre4)

    # find_min_pericentre_time()

    # plot_pericentre_info_cb()
    plot_pericentre_info()
    # plot_separation_pericentre_info()
    # plot_separation_time_info()

    plt.show()

    if len(usable_lines) > 0:
        print("\nUsable conditions in the file are on lines:", usable_lines, "\n")
    else:
        print("\nNo usable conditions were found in the simulation.\n")


if __name__ == '__main__':
    main()
