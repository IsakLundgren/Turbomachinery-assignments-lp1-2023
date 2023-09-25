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
    alpha1 = alpha3 = arctand((1 / phi) * (1 - R - psi/2))  # Eq 4.14
    alpha2 = arctand(psi/phi + tand(alpha1))
    beta2 = arctand((1 - psi) / phi - tand(alpha1))
    beta1 = arctand(psi / phi + tand(beta2))

    # Calculate pitch-chord ratios
    s_l_rotor = (DF + cosd(beta1) / cosd(beta2) - 1) / (1 / 2 * (tand(beta1) - tand(beta2)) * cosd(beta1))
    s_l_stator = (DF + cosd(alpha2) / cosd(alpha3) - 1) / (1 / 2 * (tand(alpha2) - tand(alpha3)) * cosd(alpha2))

    # Fetch profile loss coefficients
    Y_p_rotor = profileLoss(s_l_rotor, beta1, beta2)
    Y_p_stator = profileLoss(s_l_stator, alpha2, alpha3)

    # Construct total-to-total efficiency
    result = 1 - 1/2 * phi ** 2 / psi * (Y_p_rotor / (cosd(beta1) ** 2) + Y_p_stator / (cosd(alpha2) ** 2))
    return result, s_l_stator, s_l_rotor


# Return answers
eta_tt_set, dump1, dump2 = totalToTotalEff(phi_set)
eta_tt, s_l_statorList, s_l_rotorList = totalToTotalEff(phi_list)

# Print set answer
print(f'η_tt: {eta_tt_set:.4g}.')

# Plot continous efficiency
fig, ax1 = plt.subplots()
lns1 = ax1.plot(phi_list, eta_tt * 100, '-r', label='Efficiency')
ax1.set_title('Total-to-total efficiency and pitch to chord ratio')
ax1.set_ylabel('η_tt [%]')
ax1.set_xlabel('Φ [-]')
ax1.grid()
ax2 = ax1.twinx()
lns2 = ax2.plot(phi_list, s_l_statorList, '-b', label="Pitch to chord")
ax2.set_ylabel('s/l')
lns = lns1 + lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)


# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/EfficiencyAndPitchToChord.png', dpi=figureDPI)

# Show plots
plt.show()

