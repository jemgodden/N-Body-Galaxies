import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

ke = []  # List of values for kinetic energy.
pe = []  # List of values for potential energy.
RW_ke = []  # List of values for kinetic energy, exclusively for backwards interaction.
RW_pe = []  # List of values for potential energy, exclusively for backwards interaction.
total_e = []  # List of values for total energy.
RW_total_e = []  # List of values for total energy, exclusively for backwards interaction.
plot_step = []  # List of each time step, numbered.


def option_checks():
    if energy_fwds_bwds and rewind:  # Other wise will show backwards data for both lines on graph.
        print("\nPlease set rewind to 'True' in IntData.py, in order to not get the same data plotted twice.")
        exit(1)


def read_files():  # Reads the energy files made by NBody Code.
    read_energy_file()
    if energy_fwds_bwds:
        RW_ke_read()
        RW_pe_read()


def read_energy_file():  # Reads the kinetic energy files made by NBody Code.
    if rewind:
        file = open("Backwards/RewindEnergies.txt", "r")
    else:
        file = open("Forwards/Energies.txt", "r")
    for line in file:
        data = line.strip().split()
        ke.append(float(data[0]))
        ke.append(float(data[1]))
    file.close()


def RW_ke_read():  # Reads the kinetic energy files made by NBody Code.
    file = open("Backwards/RewindKE.txt", "r")
    for line in file:
        data = line.strip().split()
        RW_ke.append(float(data[0]))
    file.close()


def RW_pe_read():  # Reads the potential energy files made by NBody Code.
    file = open("Backwards/RewindPE.txt", "r")
    for line in file:
        data = line.strip().split()
        RW_pe.append(float(data[0]))
    file.close()


def find_variables():  # Finds variables used for plotting.
    find_steps()
    find_total_e()
    if energy_fwds_bwds:
        RW_find_total_e()


def find_steps():  # Creates a list of numbers that represent each time step.
    for i in range(len(ke)):
        plot_step.append(i * time_step / Gyrs)


def find_total_e():  # Calculates the total energy at each time step.
    for i in range(len(ke)):
        total_e.append(ke[i] + pe[i])


def RW_find_total_e():  # Calculates the total energy at each time step.
    for i in range(len(RW_ke)):
        RW_total_e.append(RW_ke[i] + RW_pe[i])


def manipulation():  # Manipulates the all energy from the simulation.
    ke_manip()
    pe_manip()
    total_e_manip()
    cumulative_error()


def ke_manip():  # Manipulates the kinetic energy to find the different values of errors and print them.
    avg_ke = sum(ke) / no_step
    print("\nThe average KE is", avg_ke)
    print("\nThe initial KE is", ke[0])
    print("\nThe final KE is", ke[-1])

    print("\nThe error in the final KE is", (ke[-1] - ke[0])/ke[0], "%")


def pe_manip():  # Manipulates the potential energy to find the different values of errors and print them.
    avg_pe = sum(pe) / no_step
    print("\n\nThe average PE is", avg_pe)
    print("\nThe initial PE is", pe[0])
    print("\nThe final PE is", pe[-1])

    print("\nThe error in the final PE is", (pe[-1] - pe[0])/pe[0], "%")


def total_e_manip():  # Manipulates the total energy to find the different values of errors and print them.
    avg_e = sum(total_e) / no_step
    print("\n\nThe average total energy is", avg_e)
    print("\nThe initial total energy is", total_e[0])
    print("\nThe final total energy is", total_e[-1])

    print("\nThe error in the final total energy is", (total_e[-1] - total_e[0])/total_e[0], "%")


def cumulative_error():  # Calculates the cumulative total energy error and prints it.
    c_err = 0

    E_0 = abs(ke[0]) + abs(pe[0])

    for i in range(1, len(ke), 1):
        c_err += abs(pe[i] - pe[i-1] + ke[i] - ke[i-1]) / E_0

    print("\n\nThe cumulative error for the total energy is", c_err, ".\n")


def plot():  # Plots how KE, PE and total energy change with each time step.
    plt.figure(figsize=(6, 6))
    plt.plot(plot_step, ke, 'b--', label='KE')
    plt.plot(plot_step, pe, 'r-.', label='PE')
    plt.plot(plot_step, total_e, 'k-', label='Total E')
    plt.legend(loc='upper right')
    plt.xlabel("Time (Gyrs)")
    plt.ylabel("Energy (J)")


def plot_fwbw():  # Plots how total energy changes with each time step, on an interaction run forwards and backwards.
    plt.figure(figsize=(6, 6))
    plt.plot(plot_step, total_e, 'r-', label='Forwards TE')
    plt.plot(plot_step, RW_total_e, 'b--', label='Backwards TE')
    plt.legend(loc='upper right')
    plt.xlabel("Time (Gyrs)")
    plt.ylabel("Energy (J)")


def main():

    option_checks()

    read_files()

    find_variables()

    if energy_fwds_bwds:
        print("\nNo information is available when plotting the energy forwards and backwards on the same graph.")
        plot_fwbw()
    else:
        manipulation()
        plot()

    plt.show()


if __name__ == '__main__':
    main()
