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

pericentre = []
time_of_pericentre = []
separation = []
usable_lines = []


def read_file():
    file = open("Automation_Data1.txt", "r")
    line_no = 0
    for line in file:
        line_no += 1
        if line[0] == '#':
            continue
        else:
            data = line.strip().split("\t")
            if float(data[12]) < 0.1:
                continue
            else:
                pericentre.append(float(data[12]))
                time_of_pericentre.append(-float(data[13]))
                separation.append(float(data[14]))
                if 0.2 < -float(data[13]) < 0.3 and float(data[12]) < 50:
                    usable_lines.append(line_no)
    file.close()


def plot_pericentre_info():
    plt.figure(figsize=(6, 6))
    plt.plot(pericentre, time_of_pericentre, 'k.')
    plt.gca().set_xlim([0, None])
    plt.gca().set_ylim([0, None])
    plt.fill_between([0, 50], [0.3, 0.3], [0.2, 0.2], fc='red', alpha=0.5)
    plt.xlabel("Pericentre (kpc)")
    plt.ylabel("Lookback Time (Gyrs)")


def plot_separation_pericentre_info():
    plt.figure(figsize=(6, 6))
    plt.plot(separation, pericentre, 'k.')
    plt.gca().set_ylim([0, None])
    plt.xlabel("Initial Separation (kpc)")
    plt.ylabel("Pericentre (kpc)")


def plot_separation_time_info():
    plt.figure(figsize=(6, 6))
    plt.plot(separation, time_of_pericentre, 'k.')
    plt.gca().set_ylim([0, None])
    plt.xlabel("Initial Separation (kpc)")
    plt.ylabel("Lookback Time of Pericentre (Gyrs)")


def main():  # Calling all functions in order.

    read_file()

    plot_pericentre_info()
    plot_separation_pericentre_info()
    plot_separation_time_info()

    plt.show()

    if len(usable_lines) > 0:
        print("\nUsable conditions in the file are on lines:", usable_lines)
    else:
        print("\nNo usable conditions were found in the simulation.")


if __name__ == '__main__':
    main()