# Imports
import numpy as np
import matplotlib.pyplot as plt

# Known information
alpha1p = 48  # deg
alpha2p = 16  # deg
l = 0.16  # m
i = 0  # deg
a = 0.4 * l  # m

# Set function domain
s_l = np.linspace(0.5, 2.5, 100)

# Simple angle calculations
alpha1 = alpha1p + i  # deg
camber = alpha1p - alpha2p  # deg

# Deviation at nominal incidence
m = 0.23 * (2 * a / l) ** 2 + alpha2p / 500  # -
deviation = m * camber * np.sqrt(s_l)  # deg

# Alpha2 and deflection
alpha2 = alpha2p + deviation
deflection = alpha1 - alpha2

# Plot deflection
fig, ax = plt.subplots()
ax.plot(s_l, deflection, '-r', label="Deflection")
ax.set_title('Deflection and profile loss')
ax.set_ylabel('Deflection [deg]')
ax.set_xlabel('s/l [-]')
ax.grid()
ax.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/DeflectionAndProfileLoss.png', dpi=figureDPI)
plt.show(block=True)
