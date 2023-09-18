import numpy as np


# Define degree trigonometric funtions
def sind(x):
    return np.sin(x * np.pi / 180)


def cosd(x):
    return np.cos(x * np.pi / 180)


def tand(x):
    return np.tan(x * np.pi / 180)


def arctand(x):
    return 180 / np.pi * np.arctan(x)


def profileLoss(s_l, alpha_in, alpha_out):
    # Lieblein equivalent diffusion factor
    eqDiff = (cosd(alpha_out) / cosd(alpha_in) *
              (1.12 + 0.61 * s_l * cosd(alpha_in) ** 2 * (tand(alpha_in) - tand(alpha_out))))

    # Lieblein empirical momentum thickness to coord lenght
    momThick_l = 0.004 / (1 - 1.17 * np.log(eqDiff))

    return (2 * momThick_l /
            (s_l * cosd(alpha_out) ** 2 * np.sqrt(1 / 4 * (tand(alpha_in) - tand(alpha_out)) ** 2 + 1)))


# Known quantities
phi_set = 0.7156
phi = np.linspace(0.5, 1, 100)
psi = 0.5894
DF = 0.45
Reaction = 0.5
T01 = 300  # K
P01 = 101325  # Pa
rho = 1.225  # kg m-3
gamma = 1.4
R = 287
