r"""Module to generate input velocity power spectra.

The velocity power spectrum is defined as :
$$ \frac{{\rm d} \langle v^2 \rangle}{{\rm d} \ln k} =
C \frac{k^p}{(k_{peak}^s + k^s)^{(q-p)/s}} \exp\left(-\frac{k}{k_{max}}\right) $$

Constants
=========
- RMS_VELOCITY: float
    Root mean square velocity of the fluid
- IR_SLOPE: float
    Slope of the power spectrum in the infrared, $p$ in the formula
- UV_SLOPE: float
    Slope of the power spectrum in the ultraviolet, $q$ in the formula
- SKEW: float
    Parameter $s$ in the formula
- K_MAX: float
    Wavenumber associated with the ultraviolet cutoff
- LATTICE_SIZE: int
    Number of lattice sites in one direction (assumes Lx == Ly == Lz)
- LATTICE_SPACING: float
    Size of a given lattice site
- CUTOFF: float
    Hard cutoff, just in case
"""

import numpy as np
import scipy.integrate as integrate

# Constants of the fluid
RMS_VELOCITY = 0.1
IR_SLOPE, UV_SLOPE = 5, -2 / 3
K_PEAK = 2.5e-2
SKEW = 4
K_MAX = np.pi / 4

# Constants of the mock numerical simulation
LATTICE_SIZE = 4096
LATTICE_SPACING = 0.5
CUTOFF = 2 * np.pi


def raw_power_spectrum(wave_number):
    """Raw power spectrum, yet to be normalized."""
    return (
        wave_number ** IR_SLOPE
        * (K_PEAK ** SKEW + wave_number ** SKEW) ** ((UV_SLOPE - IR_SLOPE) / SKEW)
        * np.exp(-((wave_number / K_MAX) ** 2))
    )


normalization_constant = (
    RMS_VELOCITY ** 2
    / integrate.quad(
        lambda wave_number: raw_power_spectrum(wave_number) / wave_number, 0, CUTOFF
    )[0]
)

integral_scale = (
    integrate.quad(
        lambda wave_number: raw_power_spectrum(wave_number) / wave_number ** 2,
        0,
        CUTOFF,
    )[0]
    / integrate.quad(
        lambda wave_number: raw_power_spectrum(wave_number) / wave_number, 0, CUTOFF
    )[0]
)


def normalized_power_spectrum(wave_number):
    return normalization_constant * raw_power_spectrum(wave_number)


if __name__ == "__main__":

    print("Reynolds number : {:1.1f}".format(integral_scale / LATTICE_SPACING))

    data = np.empty((LATTICE_SIZE, 3))
    for i in range(LATTICE_SIZE):
        wave_number = 2 * np.pi / LATTICE_SIZE / LATTICE_SPACING * (1 / 2 + i)
        data[i, 0] = wave_number
        data[i, 1] = normalized_power_spectrum(wave_number)
        data[i, 2] = int(np.pi * (i * 2) ** 2)

    filename = "dbpl_p-{:1.1f}_q-{:1.1f}_vrms-{:1.1e}_xi-{:1.1f}.txt".format(
        IR_SLOPE, UV_SLOPE, RMS_VELOCITY, integral_scale
    )

    np.savetxt(filename, data, fmt=["%.8f", "%.8e", "%.i"])
