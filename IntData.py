# Simulation options:
rewind = False  # Option to run the simulation backwards.
initial_txt = False  # Option to read the initial conditions from a text file.
secondary_disk = True  # Option to include a disk of particles as part of the Secondary Galaxy.
calc_energy = True  # Option to calculate energy during the simulation.
centre_mid = True  # Option to view the centre, between the two galaxies, of the interaction.

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
yg1 = 0  # y position of primary galaxy.
zg1 = 0  # y position of primary galaxy.
vxg1 = 0  # x velocity of primary galaxy.
vyg1 = 110 * km_s  # y velocity of primary galaxy.
vzg1 = 0  # z velocity of primary galaxy.
path1 = "b-"  # Colour of path for primary galaxy being plotted.

# Secondary galaxy starting conditions:
mg2 = 1e11 * sm  # Mass of secondary galaxy.
xg2 = 120 * kpc  # x position of secondary galaxy.
yg2 = 500 * kpc  # y position of secondary galaxy.
zg2 = 0  # z position of secondary galaxy.
vxg2 = 0  # x velocity of secondary galaxy.
vyg2 = -110 * km_s  # y velocity of secondary galaxy.
vzg2 = 0  # z velocity of secondary galaxy.
path2 = "r-"  # Colour of path for secondary galaxy being plotted.

# Secondary galaxy disk conditions:
no_rings = 4  # Number of rings in galaxy.
ring_rad = 4 * kpc  # Radius of innermost ring from galaxy centre.
no_rp = 4  # Number of particles in innermost ring of galaxy.
tot_rp = no_rp * (sum(range(0, no_rings + 1)))  # Total number of particles in a disk of galaxy.

# Total number of particles in simulation:
tot_part = tot_rp + 2  # Total number of particles in the simulation.
if secondary_disk:
    tot_part = tot_part + tot_rp  # Total number of particles in the simulation.

'''
To-Do list:
-Put all comments in blocks above blocks of code, unless relating to a particular line.
-Put everything into classes. Move classes out of main file, into their own one. Research classes.
-Split up functions, so each only does one thing and they all flow together.
'''
