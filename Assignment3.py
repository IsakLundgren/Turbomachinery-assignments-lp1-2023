import numpy as np
import matplotlib.pyplot as plt


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
phi_list = np.linspace(0.5, 1, 100)


def totalToTotalEff(phi):
    psi = 0.5894
    DF = 0.45
    R = 0.5

    # Calculate angles from R, psi, phi
    alpha1 = alpha3 = arctand((1 / phi) * (R + psi/2 - 1))  # Eq 4.14
    alpha2 = arctand(2 / phi * (1-R) + tand(alpha1))  # Eq 4.13a
    beta2 = arctand(tand(alpha2) - 1 / phi)  # Eq 4.15
    beta1 = arctand(2*R/phi - tand(beta2))

    # Calculate pitch-chord ratios
    s_l_rotor = (DF + cosd(alpha1) / cosd(alpha2) - 1) / (1 / 2 * (tand(alpha1) - tand(alpha2)) * cosd(alpha1))
    s_l_stator = (DF + cosd(alpha2) / cosd(alpha3) - 1) / (1 / 2 * (tand(alpha2) - tand(alpha3)) * cosd(alpha2))

    # Bugfixing
    c_x = 10
    c_1 = c_x / cosd(alpha1)
    c_2 = c_x / cosd(alpha2)
    c_3 = c_x / cosd(alpha3)
    c_t1 = c_x * tand(alpha1)
    c_t2 = c_x * tand(alpha2)
    c_t3 = c_x * tand(alpha3)

    DF_recalc_rotor = (1 - c_2 / c_1) + (c_t1 - c_t2) / (2 * c_1) * s_l_rotor
    DF_recalc_stator = (1 - c_3 / c_2) + (c_t2 - c_t3) / (2 * c_2) * s_l_stator

    # Fetch profile loss coefficients
    Y_p_rotor = profileLoss(s_l_rotor, alpha1, alpha2)
    Y_p_stator = profileLoss(s_l_stator, alpha2, alpha3)

    # Construct total-to-total efficiency
    result = 1 - 1/2 * phi ** 2 / psi * (Y_p_rotor / (cosd(beta1) ** 2) + Y_p_stator / (cosd(alpha2) ** 2))
    return result, s_l_stator, s_l_rotor


# Return answers
eta_tt_set, dump1, dump2 = totalToTotalEff(phi_set)
eta_tt, s_l_statorList, s_l_rotorList = totalToTotalEff(phi_list)

# Print set answer
print(f'η_tt: {eta_tt_set:.3g}.')

# Plot continous efficiency
fig1, ax1 = plt.subplots()
ax1.plot(phi_list, eta_tt * 100, '-r')
ax1.set_title('Total-to-total efficiency for different flow coefficients')
ax1.set_ylabel('η_tt [%]')
ax1.set_xlabel('Φ [-]')
ax1.grid()

# Plot continous rotor/stator pitch-chord ratios
fig2, ax2 = plt.subplots()
ax2.plot(phi_list, s_l_statorList, '-r', label="stator")
ax2.plot(phi_list, s_l_rotorList, '-b', label="rotor")
ax2.set_title('Pitch-chord ratios for different flow coefficients')
ax2.set_ylabel('s/l')
ax2.set_xlabel('Φ [-]')
ax2.grid()
ax2.legend()

# Show plots
plt.show()

