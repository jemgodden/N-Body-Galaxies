# Simulation options:
initial_txt = False  # Option to read the initial conditions from a text file.
galaxy_files = False  # Option to read in all galaxy particles from specific text files for each galaxy.
primary_isolation = True  # Option to run a simulation of the primary galaxy, in isolation, at (0, 0, 0).
secondary_isolation = False  # Option to run a simulation of the secondary galaxy, in isolation, at (0, 0, 0).
rewind = False  # Option to run the simulation backwards.
calc_energy = False  # Option to calculate energy during the simulation.
softening = False  # Option to include softening in the simulation.
random_disks = True

# Galaxy options:
primary_gal = True  # Option to include the primary galaxy in the interaction.
primary_disk = True  # Option to include a disk of particles as part of the Primary Galaxy.
primary_dmh_potential = True  # Option to include a dark matter halo for the primary galaxy.
primary_live_dmh = False
primary_dynamical_friction = False
secondary_gal = True  # Option to include the secondary galaxy in the interaction.
secondary_disk = True  # Option to include a disk of particles as part of the Secondary Galaxy.
secondary_dmh_potential = True  # Option to include a dark matter halo for the secondary galaxy.
secondary_live_dmh = False
secondary_dynamical_friction = False

# Simulation viewing options:
centre_mid = False  # Option to view the centre, between the two galaxies, of the interaction.
centre_pri = False  # Option to view the Primary Galaxy, during the interaction.
centre_sec = False  # Option to view the secondary Galaxy, during the interaction.
origin_pri = False  # Option to view the Primary Galaxy staying at the origin during the interaction.
origin_sec = False  # Option to view the Secondary Galaxy staying at the origin during the interaction.
gal_sep_plot = False  # Option to plot the distance between galaxies at each time step of the interaction.
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
time_step = 2e5 * yr  # Time between each step.
time_run = 0.1 * Gyr  # Total time simulation is run for.
if rewind:
    time_step = - time_step  # Time between each step, running backwards.
    time_run = - time_run  # Total time simulation is run for, backwards.
images = 11  # Total number of images shown. Remember to add 1 to make sure you get the beginning and end image.
frames = images - 1  # Number of intervals between images being shown. Could be: frames = time_run / interval.
no_step = time_run / time_step  # Total number of steps in simulation.
interval = no_step / frames  # Number of data points between image data points.
image_time_step = time_run / (frames * Gyr)  # Time between images being shown.

# Identifiers for each galaxy, used in the code to specify each.
primary = 1
secondary = 2

# Primary galaxy starting conditions:
pri_galaxy_name = "M31"  # Name of the primary galaxy.
mg1 = 2e11 * sm  # Mass of primary galaxy.
xg1 = 0 * kpc  # x position of primary galaxy.
yg1 = 0 * kpc  # y position of primary galaxy.
zg1 = 0 * kpc  # y position of primary galaxy.
vxg1 = 0 * km_s  # x velocity of primary galaxy.
vyg1 = 0 * km_s  # y velocity of primary galaxy.
vzg1 = 0 * km_s  # z velocity of primary galaxy.
pri_galaxy_marker = "bo"  # Colour and size of marker for primary galaxy being plotted.
pri_path = "b-"  # Colour of path for primary galaxy being plotted.

# Primary galaxy disk conditions:
norm_spin1 = [0, 0, 1]  # Normalised spin of the primary galaxy for x, y and z directions.
dr1 = 20 * kpc  # Radius of disk if primary galaxy.
no_rings1 = 6  # Number of rings in primary galaxy.
ring_rad1 = dr1 / no_rings1  # Radius of innermost ring from primary galaxy centre.
no_rp1 = 4  # Number of particles in innermost ring of primary galaxy.
tot_dp1 = 600  # no_rp1 * (sum(range(0, no_rings1 + 1)))  # Total number of particles in a disk of primary galaxy.
mdp1 = 7e7 * sm / tot_dp1  # Mass of each particle in the primary disk.
mu1 = 0.3  # Mean value for the gaussian distribution of particles in the primary galaxy.
sigma1 = 0.3  # Variance value for the gaussian distribution of particles in the primary galaxy.
pri_disk_name = "pDisk"  # Name of all the particles in the primary galaxy disk.
pri_disk_marker = "c."  # Colour and size of marker for primary galaxy disk particle being plotted.

# Primary galaxy dark matter halo conditions:
M_vir1 = 1.8e12 * sm  # Virial mass of the primary galaxy's dark matter halo.
R_s1 = 41 * kpc  # Scale radius of primary galaxy's dark matter halo.
c1 = 28  # Concentration of primary galaxy's dark matter halo.
R_vir1 = c1 * R_s1  # Virial radius of primary galaxy's dark matter halo.
rho_zero1 = 136 * sm / (kpc ** 3)
V_max1 = 260 * km_s
pri_dmh_name = "pDMH"

# Secondary galaxy starting conditions:
sec_galaxy_name = "M33"  # Name of the secondary galaxy.
mg2 = 8.2e10 * sm  # Mass of secondary galaxy.
xg2 = -278.1 * kpc  # x position of secondary galaxy.
yg2 = 326.7 * kpc  # y position of secondary galaxy.
zg2 = 44.5 * kpc  # z position of secondary galaxy.
vxg2 = 56.5 * km_s  # x velocity of secondary galaxy.
vyg2 = -6.5 * km_s  # y velocity of secondary galaxy.
vzg2 = 32.8 * km_s  # z velocity of secondary galaxy.
gal_marker2 = "ro"  # Colour and size of marker for primary galaxy being plotted.
path2 = "r-"  # Colour of path for secondary galaxy being plotted.

# Secondary galaxy disk conditions:
norm_spin2 = [0.048, 0.998, -0.038]  # Normalised spin of the primary galaxy for x, y and z directions.
dr2 = 9.37 * kpc  # Radius of disk if secondary galaxy.
no_rings2 = 6  # Number of rings in secondary galaxy.
ring_rad2 = dr2 / no_rings2  # Radius of innermost ring from secondary galaxy centre.
no_rp2 = 4  # Number of particles in innermost ring of secondary galaxy.
tot_dp2 = no_rp2 * (sum(range(0, no_rings2 + 1)))  # Total number of particles in a disk of secondary galaxy.
mdp2 = 1  # 7e9 * sm / tot_dp2 # Mass of each particle in the secondary disk.
mu2 = 0.3  # Mean value for the gaussian distribution of particles in the secondary galaxy.
sigma2 = 0.3  # Variance value for the gaussian distribution of particles in the secondary galaxy.
sec_disk_name = "sDisk"  # Name of all the particles in the secondary galaxy disk.
sec_disk_marker = "m."  # Colour and size of marker for secondary galaxy disk particle being plotted.

# Secondary galaxy dark matter halo conditions:
M_vir2 = 4.38e11 * sm  # Virial mass of the secondary galaxy's dark matter halo.
R_s2 = 5.0 * kpc  # Scale radius of secondary galaxy's dark matter halo.
c2 = 11  # Concentration of secondary galaxy's dark matter halo.
R_vir2 = c2 * R_s2  # Virial radius of secondary galaxy's dark matter halo.
rho_zero2 = 90 * sm / (kpc ** 3)
V_max2 = 125 * km_s
sec_dmh_name = "sDMH"

# Total number of particles in simulation:
tot_part = 0
if primary_isolation:
    tot_part += tot_dp1 + 1
if secondary_isolation:
    tot_part += tot_dp2 + 1
if not primary_isolation and not secondary_isolation:
    if primary_gal:
        tot_part += 1  # Total number of particles in the simulation.
    if primary_disk:
        tot_part += tot_dp1  # Total number of particles in the simulation.
    if secondary_gal:
        tot_part += 1  # Total number of particles in the simulation.
    if secondary_disk:
        tot_part += tot_dp2  # Total number of particles in the simulation.

# Setting the softening parameter for the simulation:
soft_param = 0.98 * (tot_part ** -0.26)

# Forcefully sets some initial conditions if running simulation of primary galaxy in isolation:
if primary_isolation:
    secondary_gal = False
    secondary_disk = False
    secondary_dmh_potential = False
    secondary_live_dmh = False
    secondary_dynamical_friction = False
    xg1 = 0
    yg1 = 0
    zg1 = 0
    vxg1 = 0
    vyg1 = 0
    vzg1 = 0

# Forcefully sets some initial conditions if running simulation of secondary galaxy in isolation:
if secondary_isolation:
    primary_gal = False
    primary_disk = False
    primary_dmh_potential = False
    primary_live_dmh = False
    primary_dynamical_friction = False
    xg2 = 0
    yg2 = 0
    zg2 = 0
    vxg2 = 0
    vyg2 = 0
    vzg2 = 0

'''
To-Do list:
-Put in Dynamical Friction.
-Put in Barnes-Hut algorithm.
-Put in a way of reading in DMH files.

-Clean up all other files.

-Put all comments in blocks above blocks of code, unless relating to a particular line.
-Put everything into classes. Move classes out of main file, into their own one. Research classes.
'''