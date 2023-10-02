import numpy as np
import statistics
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
r_m = 0.4  # m
h = 0.2  # m
m_dot = 10  # kg s-1
phi_m = 0.7156
psi_m = 0.5894
R_m = 0.5
DF_m = 0.45
h_l = 2
a_l = 0.5
rho = 1.2  # kg m-3

# Construct radial distribution vector
r = np.linspace(r_m - h/2, r_m + h/2, 1000)

# Calculate angles from R, psi, phi
alpha1_m = alpha3_m = arctand((1 / phi_m) * (1 - R_m - psi_m/2))  # Eq 4.14
alpha2_m = arctand(psi_m/phi_m + tand(alpha1_m))
beta2_m = arctand((1 - psi_m) / phi_m - tand(alpha1_m))
beta1_m = arctand(psi_m / phi_m + tand(beta2_m))

# Calculate pitch-chord ratios
s_l_rotor_m = (DF_m + cosd(beta1_m) / cosd(beta2_m) - 1) / (1 / 2 * (tand(beta1_m) - tand(beta2_m)) * cosd(beta1_m))
s_l_stator_m = (DF_m + cosd(alpha2_m) / cosd(alpha3_m) - 1) / (1 / 2 * (tand(alpha2_m) - tand(alpha3_m)) * cosd(alpha2_m))

# Fetch profile loss coefficients
Y_p_rotor_m = profileLoss(s_l_rotor_m, beta1_m, beta2_m)
Y_p_stator_m = profileLoss(s_l_stator_m, alpha2_m, alpha3_m)

# Construct total-to-total efficiency
eta_tt = 1 - 1/2 * phi_m ** 2 / psi_m * (Y_p_rotor_m / (cosd(beta1_m) ** 2) + Y_p_stator_m / (cosd(alpha2_m) ** 2))

# Calculate axial velocity, constant in radial direction
c_x = m_dot / (rho * np.pi * (r[-1] ** 2 - r[0] ** 2))  # m s-1

# Calculate K for all stations
K_1 = r_m * c_x * tand(alpha1_m)
K_2 = r_m * c_x * tand(alpha2_m)
K_3 = r_m * c_x * tand(alpha3_m)

# Calculate absolute flow angles
alpha1 = arctand(K_1/(r * c_x))
alpha2 = arctand(K_2/(r * c_x))
alpha3 = arctand(K_3/(r * c_x))

# Calculate flow angular velocity and absolute velocity
c_theta_1 = K_1 / r
c_theta_2 = K_2 / r
c_theta_3 = K_3 / r
c1 = np.sqrt(c_x ** 2 + c_theta_1 ** 2)
c2 = np.sqrt(c_x ** 2 + c_theta_2 ** 2)
c3 = np.sqrt(c_x ** 2 + c_theta_3 ** 2)

# Calculate blade velocity
U_m = c_x / phi_m
omegaRot = U_m / r_m  # Rad s-1
U = omegaRot * r  # m s-1

# Calculate relative flow angles
beta1 = arctand((U - c_theta_1)/c_x)
beta2 = arctand((U - c_theta_2)/c_x)

# Calculate blade, camber and stagger angles
alpha1_p = beta1
alpha2_p_r = beta2
alpha2_p_s = alpha2
alpha3_p = alpha3

camber_r = alpha1_p - alpha2_p_r
camber_s = alpha2_p_s - alpha3_p

stagger_r = (alpha1_p + alpha2_p_r) / 2
stagger_s = (alpha2_p_s + alpha3_p) / 2

# Calculate flow coefficient
phi = c_x / U

# Calculate stage loading
psi = phi * (tand(alpha2) - tand(alpha1))  # 5.17b

# Calculate degree of reaction
R = phi/2 * (tand(beta1) + tand(beta2))  # 5.21

# Calculate pitch to chord distirbution
l_chord = h / h_l
s_r_m = s_l_rotor_m * l_chord
s_s_m = s_l_stator_m * l_chord
tan_bha_r = s_r_m/(2 * r_m)
tan_bha_s = s_s_m/(2 * r_m)
s_r = 2 * r * tan_bha_r
s_s = 2 * r * tan_bha_s
s_l_r = s_r / l_chord
s_l_s = s_s / l_chord

# Calculate diffusion factor
DF_r = 1 - cosd(beta1) / cosd(beta2) + 1/2 * (tand(beta1) - tand(beta2)) * cosd(beta1) * s_l_r
DF_s = 1 - cosd(alpha2) / cosd(alpha3) + 1/2 * (tand(alpha2) - tand(alpha3)) * cosd(alpha2) * s_l_r

# Calculate deHaller number
deHaller_r = cosd(beta1) / cosd(beta2)
deHaller_s = cosd(alpha2) / cosd(alpha3)

# Plot angle distributions
fig, ax = plt.subplots()
ax.plot(r, alpha1, label='alpha_1')
ax.plot(r, camber_r, label='theta_rotor')
ax.plot(r, camber_s, label='theta_stator')
ax.plot(r, stagger_r, label='xi_rotor')
ax.plot(r, stagger_s, label='xi_stator')
ax.set_title('Flow, camber and stagger angles')
ax.set_ylabel('degrees [deg]')
ax.set_xlabel('r [m]')
ax.grid()
ax.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/AngleDistributions.png', dpi=figureDPI)

# Plot dimensionless quantities
fig, ax = plt.subplots()
ax.plot(r, phi, label='Flow coefficient')
ax.plot(r, psi, label='Stage load')
ax.plot(r, R, label='Degree of reaction')
ax.set_title('Normal stage parameters')
ax.set_ylabel('[-]')
ax.set_xlabel('r [m]')
ax.grid()
ax.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/NormalStageParameters.png', dpi=figureDPI)

# Plot deHaller and diffusion factor
fig, ax = plt.subplots()
ax.plot(r, DF_r, label='Diffusion factor rotor')
ax.plot(r, DF_s, label='Diffusion factor stator')
ax.plot(r, deHaller_r, label='DeHaller number rotor')
ax.plot(r, deHaller_s, label='DeHaller number stator')
ax.set_title('Diffusion factor and DeHaller number')
ax.set_ylabel('[-]')
ax.set_xlabel('r [m]')
ax.grid()
ax.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/DFandDeHaller.png', dpi=figureDPI)

# Show plots
plt.show()

