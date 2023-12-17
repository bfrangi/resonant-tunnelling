from numpy import seterr
seterr(all='raise')


####################################################################################################
#                                        PotentialSimulator                                        #
####################################################################################################


class PotentialSimulator:
    def __init__(self, V: callable):
        self.V = V
        self.x_max = 2.0  # nm
        self.x_density = 100
        self.x_arr = None
        self.V_arr = None
        self.plt = None

    def compute_spatial_array(self):
        from numpy import array
        x_max = self.x_max
        density = self.x_density
        self.x_arr = array([i/density for i in range(int(density*x_max)+1)])
        return self.x_arr

    def compute_potential_array(self):
        from numpy import array
        x_arr = self.x_arr if self.x_arr is not None else self.compute_spatial_array()
        self.V_arr = array([self.V(x) for x in x_arr])
        return self.V_arr

    def generate_potential_plot(self):
        import matplotlib.pyplot as plt
        V_arr = self.V_arr if self.V_arr is not None else self.compute_potential_array()
        x_arr = self.x_arr
        plt.plot(x_arr, V_arr)
        plt.xlim([x_arr[0], x_arr[-1]])
        plt.ylim([min(V_arr), 1.2*max(V_arr)])
        plt.xlabel("x (nm)")
        plt.ylabel("V (eV)")
        plt.title("Potential")
        self.plt = plt
        return self.plt

    def show_potential_plot(self):
        plt = self.plt if self.plt is not None else self.generate_potential_plot()
        plt.show()


####################################################################################################
#                                       TransmissionSimulator                                      #
####################################################################################################


class TransmissionSimulator:
    def __init__(self, potential_simulator: PotentialSimulator):
        self.potential_simulator = potential_simulator
        self.E_density = 1000
        self.E_top_margin = 0.5
        self.E_arr = None
        self.T_arr = None
        self.x_arr = None
        self.V_arr = None
        self.plt = None

        from numpy import sqrt

        # Unit Base
        self.m_e = 9.1093837e-31  # kg = M
        self.eV = 1.602176634e-19  # J = ML²T⁻²
        self.nm = 1e-9  # m = L

        # Adimensionalised constants
        hbar = 1.0545718e-34  # J·s = ML²T⁻¹
        self.hbar = hbar/(sqrt(self.eV)*self.nm*sqrt(self.m_e)) # eV^(1/2) nm m_e^(1/2)
        self.m = 1.0  # m_e

    def compute_energy_array(self):
        from numpy import array
        V_arr = self.potential_simulator.compute_potential_array()
        max_E = max(V_arr) * (1 + self.E_top_margin)
        density = self.E_density
        max_i = int(max_E * density)
        self.E_arr = array([i/density for i in range(1, max_i)])
        return self.E_arr

    # Remove points with k = 0
    @staticmethod
    def clean_arrays(x_arr, k_arr, verbose=False):
        from numpy import delete
        indices = []
        for i in range(len(k_arr)):
            if not k_arr[i]:
                indices.append(i)
        x_arr = delete(x_arr, indices)
        k_arr = delete(k_arr, indices)
        if verbose:
            if indices:
                print(f"Removed some points:", indices)
            else:
                print(f"Nothing removed")
        return x_arr, k_arr

    def compute_transmission(self, E, verbose=False):
        from numpy import array, absolute, exp, emath
        if self.V_arr is None:
            self.V_arr = self.potential_simulator.compute_potential_array()
        V_arr = self.V_arr
        if self.x_arr is None:
            self.x_arr = self.potential_simulator.compute_spatial_array()
        x_arr = self.x_arr

        # Initial propagation matrix
        P = array([
            [1, 0],
            [0, 1],
        ], dtype=complex)

        # Wave number
        k_arr = array([
            emath.sqrt(2*self.m*(E-V_val))/self.hbar
            for V_val in V_arr
        ])

        # Clean arrays
        x_arr, k_arr = self.clean_arrays(x_arr, k_arr, verbose=verbose)

        # Compute propagation matrix
        for i in range(len(x_arr) - 1):
            dL = x_arr[i+1] - x_arr[i]
            frac = k_arr[i + 1]/k_arr[i]
            exponent = 1j*k_arr[i]*dL
            P @= array([
                [(1 + frac) * exp(-exponent),
                 (1 - frac) * exp(-exponent)],
                [(1 - frac) * exp(exponent),
                 (1 + frac) * exp(exponent)],
            ]) / 2

        # Compute transmission coefficient
        T = 1/absolute(P[0, 0])**2
        return T

    def compute_transmission_array(self, verbose=False):
        from numpy import array
        E_arr = self.E_arr if self.E_arr is not None else self.compute_energy_array()
        T_arr = []
        for E in E_arr:
            T_arr.append(self.compute_transmission(E, verbose=verbose))
        self.T_arr = array(T_arr)
        return self.T_arr

    def generate_transmission_plot(self):
        import matplotlib.pyplot as plt
        E_arr = self.E_arr if self.E_arr is not None else self.compute_energy_array()
        T_arr = self.T_arr if self.T_arr is not None else self.compute_transmission_array()
        plt.plot(E_arr, T_arr)
        plt.xlim([E_arr[0], E_arr[-1]])
        plt.xlabel("E (eV)")
        plt.ylabel("T")
        plt.title("Transmission coefficient as a function of the electron energy")
        self.plt = plt
        return self.plt

    def show_transmission_plot(self):
        plt = self.plt if self.plt is not None else self.generate_transmission_plot()
        plt.show()


####################################################################################################
#                       AnalyticDoubleSymmetricPotentialTransmissionSimulator                      #
####################################################################################################


class AnalyticDoubleSymmetricPotentialTransmissionSimulator:
    def __init__(self, V1=1.0, V2=1.0):
        self.V1 = V1
        self.V2 = V2
        self.w_1 = 0.2  # nm
        self.w_2 = 0.4  # nm
        self.w_3 = 0.6  # nm
        self.w_4 = 0.4  # nm
        self.w_5 = 0.4  # nm
        self.E_density = 1000
        self.E_top_margin = 0.5
        self.E_arr = None
        self.E_max = None
        self.T_arr = None
        self.plt = None

        from numpy import sqrt

        # Unit Base
        self.m_e = 9.1093837e-31  # kg = M
        self.eV = 1.602176634e-19  # J = ML²T⁻²
        self.nm = 1e-9  # m = L

        # Adimensionalised constants
        hbar = 1.0545718e-34  # J·s = ML²T⁻¹
        self.hbar = hbar/(sqrt(self.eV)*self.nm*sqrt(self.m_e)) # eV^(1/2) nm m_e^(1/2)
        self.m = 1.0  # m_e

    def compute_energy_array(self):
        from numpy import array
        max_E = self.E_max if self.E_max is not None else max(self.V1, self.V2)
        max_E = max_E * (1 + self.E_top_margin)
        density = self.E_density
        max_i = int(max_E * density)
        self.E_arr = array([i/density for i in range(1, max_i)])
        return self.E_arr
    
    def k(self, E, V):
        from numpy import emath
        return emath.sqrt(2*self.m*(E-V)/self.hbar**2)


    def X(self, E, V):
        from numpy import emath
        return emath.sqrt(2*self.m*(V-E)/self.hbar**2)


    def compute_transmission(self, E):
        from numpy import arctan, exp, absolute
        try:
            k_1 = self.k(E, 0)
            X_2 = self.X(E, self.V1)
            k_3 = self.k(E, 0)
            X_4 = self.X(E, self.V2)
            k_5 = self.k(E, 0)

            phi_1 = k_3*self.w_3
            phi_2 = arctan(X_2/k_1)
            phi_3 = arctan(X_2/k_3)
            phi_4 = arctan(X_4/k_3)
            phi_5 = arctan(X_4/k_5)

            K = exp(X_2*self.w_2 + X_4*self.w_4)*(
                exp(1j*(- phi_1 + phi_2 + phi_3 + phi_4 + phi_5))
                - exp(1j*(phi_1 + phi_2 - phi_3 - phi_4 + phi_5))
            ) + exp(X_2*self.w_2 - X_4*self.w_4)*(
                - exp(1j*(- phi_1 + phi_2 + phi_3 - phi_4 - phi_5))
                + exp(1j*(phi_1 + phi_2 - phi_3 + phi_4 - phi_5))
            ) + exp(-X_2*self.w_2 + X_4*self.w_4)*(
                - exp(1j*(- phi_1 - phi_2 - phi_3 + phi_4 + phi_5))
                + exp(1j*(phi_1 - phi_2 + phi_3 - phi_4 + phi_5))
            ) + exp(-X_2*self.w_2 - X_4*self.w_4)*(
                exp(1j*(- phi_1 - phi_2 - phi_3 - phi_4 - phi_5))
                - exp(1j*(phi_1 - phi_2 + phi_3 + phi_4 - phi_5))
            )

            return (2.0**8 * k_1 * (X_2 * k_3 * X_4)**2 * k_5) / absolute(K)**2 \
                / (k_1**2 + X_2**2) / (k_3**2 + X_2**2) / (k_3**2 + X_4**2) / (k_5**2 + X_4**2)
        except:
            return 1


    def compute_transmission_array(self):
        from numpy import array, absolute
        E_arr = self.E_arr if self.E_arr is not None else self.compute_energy_array()
        T_arr = []
        for E in E_arr:
            T_arr.append(absolute(self.compute_transmission(E)))
        self.T_arr = array(T_arr)
        return self.T_arr

    def generate_transmission_plot(self):
        import matplotlib.pyplot as plt
        E_arr = self.E_arr if self.E_arr is not None else self.compute_energy_array()
        T_arr = self.T_arr if self.T_arr is not None else self.compute_transmission_array()
        plt.plot(E_arr, T_arr)
        plt.xlim([E_arr[0], E_arr[-1]])
        plt.xlabel("E (eV)")
        plt.ylabel("T")
        plt.title("Transmission coefficient as a function of the electron energy")
        self.plt = plt
        return self.plt

    def show_transmission_plot(self):
        plt = self.plt if self.plt is not None else self.generate_transmission_plot()
        plt.show()
