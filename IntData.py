# Simulation options:
initial_txt = False  # Option to read the initial conditions from a text file.
galaxy_files = False  # Option to read in all galaxy particles from specific text files for each galaxy.
primary_isolation = False  # Option to run a simulation of the primary galaxy, in isolation, at (0, 0, 0).
secondary_isolation = False  # Option to run a simulation of the secondary galaxy, in isolation, at (0, 0, 0).
rewind = False  # Option to run the simulation backwards.
calc_energy = False  # Option to calculate energy during the simulation.

# Galaxy options:
primary_gal = True  # Option to include the primary galaxy in the interaction.
primary_disk = True  # Option to include a disk of particles as part of the Primary Galaxy.
primary_dmh = False  # Option to include a dark matter halo for the primary galaxy.
secondary_gal = False  # Option to include the secondary galaxy in the interaction.
secondary_disk = False  # Option to include a disk of particles as part of the Secondary Galaxy.
secondary_dmh = False  # Option to include a dark matter halo for the secondary galaxy.

# Simulation viewing options:
centre_mid = False  # Option to view the centre, between the two galaxies, of the interaction.
centre_pri = True  # Option to view the Primary Galaxy, during the interaction.
energy_fwds_bwds = False  # Option to view forwards and backwards energy of an interaction, on one graph.

# Constants:
G = 6.67e-11  # Gravitational constant.
km_s = 1e3  # Kilometres per second.
pc = 3.086e16  # Parsec in metres.
kpc = 1e3 * pc  # Kilo-parsec in metres.
sm = 1.989e30  # Solar Mass.
yr = 60 * 60 * 24 * 365  # Year in seconds.
Gyr = 1e9 * yr  # Giga-year in seconds.

# Simulation time conditions:
time_step = 5e6 * yr  # Time between each step.
time_run = 3 * Gyr  # Total time simulation is run for.
if rewind:
    time_step = - time_step  # Time between each step, running backwards.
    time_run = - time_run  # Total time simulation is run for, backwards.
images = 16  # Total number of images shown. Remember to add 1 to make sure you get the beginning and end image.
frames = images - 1  # Number of intervals between images being shown. Could be: frames = time_run / interval.
no_step = time_run / time_step  # Total number of steps in simulation.
interval = no_step / frames  # Number of data points between image data points.
image_time_step = time_run / (frames * Gyr)  # Time between images being shown.

# Primary galaxy starting conditions:
mg1 = 1e11 * sm  # Mass of primary galaxy.
xg1 = 0  # x position of primary galaxy.
yg1 = 0 * kpc  # y position of primary galaxy.
zg1 = 0  # y position of primary galaxy.
vxg1 = 0  # x velocity of primary galaxy.
vyg1 = 0 * km_s  # 110 * kms # y velocity of primary galaxy.
vzg1 = 0  # z velocity of primary galaxy.
path1 = "b-"  # Colour of path for primary galaxy being plotted.

# Primary galaxy disk conditions:
norm_spin1 = [0.5, 0.25, 0.829]  # Normalised spin of the primary galaxy for x, y and z directions.
dr1 = 20 * kpc  # Radius of disk if primary galaxy.
no_rings1 = 6  # Number of rings in primary galaxy.
ring_rad1 = dr1 / no_rings1  # Radius of innermost ring from primary galaxy centre.
no_rp1 = 4  # Number of particles in innermost ring of primary galaxy.
tot_rp1 = no_rp1 * (sum(range(0, no_rings1 + 1)))  # Total number of particles in a disk of primary galaxy.

# Primary galaxy dark matter halo conditions:
M_vir1 = 9e11 * sm  # Virial mass of the primary galaxy's dark matter halo.
R_s1 = 15 * kpc  # Scale radius of primary galaxy's dark matter halo.
c1 = 12  # Concentration of primary galaxy's dark matter halo.
R_vir1 = c1 * R_s1  # Virial radius of primary galaxy's dark matter halo.

# Secondary galaxy starting conditions:
mg2 = 1e11 * sm  # Mass of secondary galaxy.
xg2 = 0  # 120 * kpc  # x position of secondary galaxy.
yg2 = 0 * kpc  # 500 * kpc  # y position of secondary galaxy.
zg2 = 0  # z position of secondary galaxy.
vxg2 = 0  # x velocity of secondary galaxy.
vyg2 = 0 * km_s  # -75 * kms  # -110 * km_s  # y velocity of secondary galaxy.
vzg2 = 0  # z velocity of secondary galaxy.
path2 = "r-"  # Colour of path for secondary galaxy being plotted.

# Secondary galaxy disk conditions:
norm_spin2 = [0, 0, 1]  # Normalised spin of the primary galaxy for x, y and z directions.
dr2 = 20 * kpc  # Radius of disk if secondary galaxy.
no_rings2 = 6  # Number of rings in secondary galaxy.
ring_rad2 = dr2 / no_rings2  # Radius of innermost ring from secondary galaxy centre.
no_rp2 = 4  # Number of particles in innermost ring of secondary galaxy.
tot_rp2 = no_rp2 * (sum(range(0, no_rings2 + 1)))  # Total number of particles in a disk of secondary galaxy.

# Secondary galaxy dark matter halo conditions:
M_vir2 = 9e11 * sm  # Virial mass of the secondary galaxy's dark matter halo.
R_s2 = 15 * kpc  # Scale radius of secondary galaxy's dark matter halo.
c2 = 12  # Concentration of secondary galaxy's dark matter halo.
R_vir2 = c2 * R_s2  # Virial radius of secondary galaxy's dark matter halo.

# Total number of particles in simulation:
tot_part = tot_rp1 + 1  # Total number of particles in the simulation.
if secondary_gal:
    tot_part = tot_part + 1  # Total number of particles in the simulation.
if secondary_disk:
    tot_part = tot_part + tot_rp2  # Total number of particles in the simulation.

'''
To-Do list:
-Put all comments in blocks above blocks of code, unless relating to a particular line.
-Get rid of any repeated code, especially with x, y and z of a coordinate being used (i.e in being printed.)
-Split up functions, so each only does one thing and they all flow together.
-Put everything into classes. Move classes out of main file, into their own one. Research classes.
'''