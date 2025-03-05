from typing_extensions import Self
import math

from src.constants import *


class Body:
    def __init__(
            self,
            name: str,
            mass: float,
            position: list[float],
            velocity: list[float],
            colour: str
    ):
        self._name: str = name
        self._mass: float = mass
        self._cur_position: list[float] = position  # Array of x, y and z position of body.
        self._cur_velocity: list[float] = velocity  # Array of x, y and z velocities of body.
        self._cur_acceleration: list[float] = [0, 0, 0]  # Array of x, y and z acceleration of body.
        self._position_history: list[list[float]] = [[], [], []]  # Array of all x, y and z position values of body.
        self._velocity_history: list[list[float]] = [[], [], []]  # Array of all x, y and z velocity values of body.
        self._colour: str = colour

    @property
    def name(self) -> str:
        return self.name

    @property
    def mass(self) -> float:
        return self._mass

    @property
    def cur_position(self) -> list[float]:
        return self._cur_position

    @property
    def cur_velocity(self) -> list[float]:
        return self._cur_velocity

    @property
    def cur_acceleration(self) -> list[float]:
        return self._cur_acceleration

    @property
    def position_history(self) -> list[list[float]]:
        return self._position_history

    @property
    def velocity_history(self) -> list[list[float]]:
        return self._velocity_history

    @property
    def colour(self) -> str:
        return self._colour

    def find_separation(self, other: Self) -> tuple[list[float], float]:
        rx = self.cur_position[0] - other.cur_position[0]
        ry = self.cur_position[1] - other.cur_position[1]
        rz = self.cur_position[2] - other.cur_position[2]
        r = math.sqrt((rx ** 2) + (ry ** 2) + (rz ** 2))
        return [rx, ry, rz], r

    def find_relative_velocity(self, other: Self) -> tuple[list[float], float]:
        vx = self.cur_velocity[0] - other.cur_velocity[0]
        vy = self.cur_velocity[1] - other.cur_velocity[1]
        vz = self.cur_velocity[2] - other.cur_velocity[2]
        v = math.sqrt((vx ** 2) + (vy ** 2) + (vz ** 2))
        return [vx, vy, vz], v

    def calculate_newtonian_acceleration(self, other: Self) -> None:
        r_xyz, r = self.find_separation(other)

        if softening:
            a = - (G * other.m) / ((r ** 2) + (soft_param ** 2))
        else:
            a = - (G * other.m) / (r ** 2)

        for i in range(len(r_xyz)):
            self.cur_acceleration[i] += a * (r_xyz[i] / r)

    def calculate_dmh_acceleration(self, other: Self, galaxy_id) -> None:
        a = 0
        r_xyz, r = self.find_separation(other)

        if galaxy_id == primary:
            a = ((G * M_vir1) / (math.log(1 + c1) - (c1 / (1 + c1)))) * \
                (((r / (r + R_s1)) - (math.log(1 + r / R_s1))) / (r ** 2))
        elif galaxy_id == secondary:
            a = ((G * M_vir2) / (math.log(1 + c2) - (c2 / (1 + c2)))) * \
                (((r / (r + R_s2)) - (math.log(1 + r / R_s2))) / (r ** 2))

        for i in range(len(r_xyz)):
            self.cur_acceleration[i] += a * (r_xyz[i] / r)

    def calculate_dynamical_friction(self, other: Self, causing_galaxy_id) -> None:
        r_xyz, r = self.find_separation(other)
        v_xyz, v = self.find_relative_velocity(other)

        density_distribution = 0
        v_dispersion = 0

        # epsilon = (0.98 * (3002 ** -0.26)) * kpc
        epsilon = 28.5 * KPC
        ln_lambda = math.log(r / (1.4 * epsilon))

        if causing_galaxy_id == primary:
            density_distribution = (102 * critical_density) / ((r / R_s1) * (1 + (r / R_s1)) ** 2)
            v_dispersion = V_max1 * ((1.4393 * (r / R_s1) ** 0.354) / (1 + 1.1756 * (r / R_s1) ** 0.725))

        elif causing_galaxy_id == secondary:
            density_distribution = rho_zero2 / ((r / R_s2) * (1 + (r / R_s2)) ** 2)
            v_dispersion = V_max2 * ((1.4393 * (r / R_s2) ** 0.354) / (1 + 1.1756 * (r / R_s2) ** 0.725))

        X = v / ((2 ** 0.5) * v_dispersion)

        a = - ((4 * math.pi * (G ** 2) * self.mass * ln_lambda * density_distribution) / (v ** 2)) * (
                math.erf(X) - (2 * X / (math.pi ** 0.5)) * math.exp(-(X ** 2)))

        for i in range(len(self.cur_acceleration)):
            self.cur_acceleration[i] += a * (v_xyz[i] / v)

    def calculate_kinetic_energy(self) -> float:
        v = (self.cur_velocity[0] ** 2) + (self.cur_velocity[1] ** 2) + (self.cur_velocity[2] ** 2)

        return 0.5 * self.mass * (v ** 2)

    def calculate_potential_energy(self, other: Self) -> float:
        r_xyz, r = self.find_separation(other)

        return - (G * self.mass * other.m) / (2 * r)

    def calculate_dmh_potential_energy(self, galaxy_id) -> float:
        halo_pe = 0

        if galaxy_id == primary:
            halo_pe = -(G * M_vir1 / r) * (1 / (math.log(1 + c1) - (c1 / (1 + c1)))) * math.log(1 + (r / R_s1))
        elif galaxy_id == secondary:
            halo_pe = -(G * M_vir2 / r) * (1 / (math.log(1 + c2) - (c2 / (1 + c2)))) * math.log(1 + (r / R_s2))

        return halo_pe * (self.mass / 2)
