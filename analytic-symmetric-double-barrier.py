from numpy import sqrt, arctan, exp, absolute
import matplotlib.pyplot as plt

# Unit Base
m_e = 9.1093837e-31  # kg = M
eV = 1.602176634e-19  # J = ML²T⁻²
nm = 1e-9  # m = L

# Adimensionalised constants
hbar = 1.0545718e-34  # J·s = ML²T⁻¹
hbar /= (sqrt(eV)*nm*sqrt(m_e))  # eV^(1/2) nm m_e^(1/2)
m = 1.0  # m_e

# Regions of the potential
w_1 = 0.2  # nm
w_2 = 0.4  # nm
w_3 = 0.6  # nm
w_4 = 0.4  # nm
w_5 = 0.4  # nm

V_2 = 1.0  # eV
V_4 = 1.0  # eV

# Functions


def k(m, E, V):
    return sqrt(2*m*(E-V)/hbar**2)


def X(m, E, V):
    return sqrt(2*m*(V-E)/hbar**2)


def T(E):
    k_1 = k(m, E, 0)
    X_2 = X(m, E, V_2)
    k_3 = k(m, E, 0)
    X_4 = X(m, E, V_4)
    k_5 = k(m, E, 0)

    phi_1 = k_3*w_3
    phi_2 = arctan(X_2/k_1)
    phi_3 = arctan(X_2/k_3)
    phi_4 = arctan(X_4/k_3)
    phi_5 = arctan(X_4/k_5)

    K = exp(X_2*w_2 + X_4*w_4)*(
        exp(1j*(- phi_1 + phi_2 + phi_3 + phi_4 + phi_5))
        - exp(1j*(phi_1 + phi_2 - phi_3 - phi_4 + phi_5))
    ) + exp(X_2*w_2 - X_4*w_4)*(
        - exp(1j*(- phi_1 + phi_2 + phi_3 - phi_4 - phi_5))
        + exp(1j*(phi_1 + phi_2 - phi_3 + phi_4 - phi_5))
    ) + exp(-X_2*w_2 + X_4*w_4)*(
        - exp(1j*(- phi_1 - phi_2 - phi_3 + phi_4 + phi_5))
        + exp(1j*(phi_1 - phi_2 + phi_3 - phi_4 + phi_5))
    ) + exp(-X_2*w_2 - X_4*w_4)*(
        exp(1j*(- phi_1 - phi_2 - phi_3 - phi_4 - phi_5))
        - exp(1j*(phi_1 - phi_2 + phi_3 + phi_4 - phi_5))
    )

    return (2.0**8 * k_1 * (X_2 * k_3 * X_4)**2 * k_5) / absolute(K)**2 \
        / (k_1**2 + X_2**2) / (k_3**2 + X_2**2) / (k_3**2 + X_4**2) / (k_5**2 + X_4**2)


max_E = min(V_2, V_4)
density = 1000
max_i = int(max_E * density)
E_list = [i/density for i in range(1, max_i)]
T_list = [T(E) for E in E_list]

# Plot
plt.plot(E_list, T_list)
plt.xlim([0, E_list[-1]])
plt.ylim([0, 1])
plt.show()
