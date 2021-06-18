r"""Module to generate input velocity power spectra.

The velocity power spectrum is defined as :
$$ \frac{{\rm d} \langle v^2 \rangle}{{\rm d} \ln k} =
C \frac{k^p}{(k_{peak}^s + k^s)^{(q-p)/s}} \exp\left(-\frac{k}{k_{max}}\right) $$

Constants
=========
- RMS_VELOCITY: float
    Root mean square velocity of the fluid
- INTEGRAL_SCALE: float
    Integral scale of the fluid, in physical units
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
import scipy.optimize as opt

# Constants of the fluid
RMS_VELOCITY = 0.1
INTEGRAL_SCALE = 10
IR_SLOPE, UV_SLOPE = 5, -2 / 3
SKEW = 2
K_MAX = np.pi / 4

# Constants of the mock numerical simulation
LATTICE_SIZE = 4096
LATTICE_SPACING = 0.5
CUTOFF = 2 * np.pi


def raw_power_spectrum(wave_number, k_peak):
    """Raw power spectrum, yet to be normalized."""
    return (
        wave_number ** IR_SLOPE
        * (k_peak ** SKEW + wave_number ** SKEW)
        ** ((UV_SLOPE - IR_SLOPE) / SKEW)
        * np.exp(-((wave_number / K_MAX) ** 2))
    )

def normalized_power_spectrum(wave_number,normalization,k_peak):
    return normalization * raw_power_spectrum(
        wave_number, k_peak
    )


def square_velocity(dummy_norm, dummy_peak):

    return integrate.quad(
        lambda wave_number: normalized_power_spectrum(
            wave_number, dummy_norm, dummy_peak
        )
        / wave_number,
        0,
        CUTOFF,
    )[0]


def integral_scale(dummy_norm, dummy_peak):

    return (
        np.pi
        / 4
        * (
            integrate.quad(
                lambda wave_number: normalized_power_spectrum(
                    wave_number, dummy_norm, dummy_peak
                )
                / wave_number ** 2,
                0,
                CUTOFF,
            )[0]
            / integrate.quad(
                lambda wave_number: normalized_power_spectrum(
                    wave_number, dummy_norm, dummy_peak
                )
                / wave_number,
                0,
                CUTOFF,
            )[0]
        )
    )


normalization, k_peak = opt.root(
    lambda x: [
        INTEGRAL_SCALE - integral_scale(*x),
        RMS_VELOCITY ** 2 - square_velocity(*x),
    ],
    [1, 1],
).x




if __name__ == "__main__":

    print("Reynolds number : {:1.1f}".format(INTEGRAL_SCALE / LATTICE_SPACING))
    
    print("vrms : {:.2g} integral scale : {:.2g}".format(square_velocity(normalization, k_peak),
                                                          integral_scale(normalization, k_peak)))
    data = np.empty((LATTICE_SIZE, 3))
    for i in range(LATTICE_SIZE):
        wave_number = 2 * np.pi / LATTICE_SIZE / LATTICE_SPACING * (1 / 2 + i)
        data[i, 0] = wave_number
        data[i, 1] = normalized_power_spectrum(wave_number, normalization, k_peak)
        data[i, 2] = int(np.pi * (i * 2) ** 2)

    filename = "dbpl_p-{:1.1f}_q-{:1.1f}_vrms-{:1.1e}_xi-{:1.1f}.txt".format(
        IR_SLOPE, UV_SLOPE, RMS_VELOCITY, INTEGRAL_SCALE
    )

    np.savetxt(filename, data, fmt=["%.8f", "%.8e", "%.i"])
