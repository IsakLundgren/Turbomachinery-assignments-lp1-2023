# Imports
import numpy as np
import matplotlib.pyplot as plt

# Known information
alpha1p = 48  # rad
alpha2p = 16  # deg
l = 0.16  # m
i = 0  # deg
a = 0.4 * l  # m

# Set function domain
s_l = np.linspace(0.5, 2.75, 1000)


# Define degree trigonometric funtions
def sind(x):
    return np.sin(x * np.pi / 180)


def cosd(x):
    return np.cos(x * np.pi / 180)


def tand(x):
    return np.tan(x * np.pi / 180)


# Simple angle calculations
alpha1 = alpha1p + i  # deg
camber = alpha1p - alpha2p  # deg

# Deviation at nominal incidence
deviation = (0.23 * (2 * a / l) ** 2 + alpha2p / 500) / (1-(camber * np.sqrt(s_l)) / 500) * camber * np.sqrt(s_l)  # deg

# Alpha2 and deflection
alpha2 = alpha2p + deviation
alpha2sergio = 500 * (0.23*(2 * a/l) ** 2 * camber * s_l ** 0.5 + alpha2p) / (500 - camber * s_l ** 0.5)
deflection = alpha1 - alpha2
deflection_sergio = alpha1 - alpha2

# Plot deflection
fig, ax1 = plt.subplots()
ax1.plot(s_l, deflection, '-r', label="Deflection")
ax1.plot(s_l, deflection_sergio, '-g', label="Deflection - sergio")
ax1.set_title('Deflection and profile loss')
ax1.set_ylabel('Deflection [deg]')
ax1.set_xlabel('s/l [-]')
ax1.grid()
ax1.legend()

# Lieblein equivalent diffusion factor
eqDiff = (cosd(alpha2) / cosd(alpha1) *
          (1.12 + 0.61 * s_l * cosd(alpha1) ** 2 * (tand(alpha1) - tand(alpha2))))

# Lieblein empirical momentum thickness to coord lenght
momThick_l = 0.004 / (1 - 1.17 * np.log(eqDiff))

# Pressure loss coefficient
pressureLoss = (2 * momThick_l /
                (s_l * l * cosd(alpha2) ** 2 * np.sqrt(1/4 * (tand(alpha1) - tand(alpha2)) ** 2 + 1)))

# Plot pressure loss coefficient
ax2 = ax1.twinx()
ax2.plot(s_l, pressureLoss, '-b', label='Pressure loss coefficient')
ax2.set_ylabel('Pressure loss coefficient [-]')
ax2.legend(loc=3)

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/DeflectionAndProfileLoss.png', dpi=figureDPI)

fig = plt.figure()
# plt.plot(s_l, alpha2, '-b')
# plt.plot(s_l, alpha1*np.ones(len(s_l)), '-r')
plt.plot(s_l, momThick_l, '-g')
# Show figures
plt.show(block=True)

# Find eqDiff at s_l = 1 and s_l = 2
s_lClosestTo1 = min(list(s_l), key=lambda x: abs(1 - x))
s_lClosestTo2 = min(list(s_l), key=lambda x: abs(2 - x))
ind1 = list(s_l).index(s_lClosestTo1)
ind2 = list(s_l).index(s_lClosestTo2)

print(f'c_max by c_2 at s/l = 1: {eqDiff[ind1]:.3f}.')
print(f'c_max by c_2 at s/l = 2: {eqDiff[ind2]:.3f}.')
